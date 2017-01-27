function [O U] = lap2d2p(iprec,ns,s,ich,ch,idip,dst,dv,ipot,igr,ihe,...
                         nt,t,ipott,igrt,ihet,e1,e2,o)
% draft of doubly periodic in 2D Laplace wrapper around CMCL FMMs.
%
% Simplifying assumptions for now:
%  * no self-evals just extra targs.
%  * source strengths compatible w/ periodizing.
%  * all sources and targs lie in the unit cell.
%  * unit cell centered on origin.
%
% kernel is usual (1/2pi) log (1/r). Note this is -1/2pi times the CMCL
% normalization. Coordinates are as stacks of 2-cmpt real col vecs.
% Sources, dipole strengths are real valued.
%
% unit cell centered at 0 with lattice vecs e1 and e2 each expressed as complex
% numbers. Last arg o is optional options struct.
% All other args, and output struct O, are as for lfmm2dpart.
%
% Barnett, afternoon of 1/26/17
if nargin==0, test_lap2d2p; return; end
if nargin<19, o=[]; end
if ~isfield(o,'verb'), o.verb=0; end     % how much diagnostic output
if ~isfield(o,'noperi'), o.noperi=0; end
badness = cond([real([e1 e2]);imag([e1 e2])]); % cond(lattice vector matrix)
if badness>10, warning('unit cell close to singular, may fail!'); end
cmcl = -2*pi;       % factor by which CMCL normalization is off
ch = ch(:)'*(1/cmcl); dst = dst(:)'*(1/cmcl);   % correct & ensure row vecs

% walls
U.e1 = e1; U.e2 = e2;
m = 25; [U L R B T] = doublywalls(U,m);  % note U is unit cell not CMCL U.
n1 = T.nx(1); n2 = R.nx(1);           % the two wall normals
% proxy pts
M = ceil(80*badness);          % # proxy
Rp = 1.4; p.c = Rp * exp(1i*(0:M-1)'/M*2*pi);   % circle
p.x = e1*real(p.c) + e2*imag(p.c);     % ellipse
p = setupquad(p);         % adds s.w weights and s.nx normals to proxy
if o.verb>1, figure; showsegment({L R T B}); hold on; plot(p.x,'r.'); end
% periodizing Q block
Q = Qmat(p,L,R,B,T,@LapSLP);    % use proxy monopoles

% compute discrep: one src to 12 targ wall copies (ns -> around 300 targs)
tw = [R.x+e1-e2;R.x+e1;R.x+e1+e2;L.x-e1-e2;L.x-e1;L.x-e1+e2;
     T.x-e1+e2;T.x+e2;T.x+e1+e2;B.x-e1-e2;B.x-e2;B.x+e1-e2].';  % row vec in C
tw = [real(tw);imag(tw)];   % 2-by-n format
W = lfmm2dpart(iprec,ns,s,ich,ch,idip,dst,dv,0,0,0,12*m,tw,1,1,0);
f = discrepblk(W.pottarg(1:6*m));       % 1st discrep block
fn = discrepblk([real(n2),imag(n2)] * W.gradtarg(:,1:6*m));  % wall n-derivs
g = discrepblk(W.pottarg(6*m+1:end));
gn = discrepblk([real(n1),imag(n1)] * W.gradtarg(:,6*m+1:end));
d = -real([f,fn,g,gn])';     % cancel out discrep col vec

warning('off','MATLAB:nearlySingularMatrix')  % backward-stable ill-cond is ok!
warning('off','MATLAB:rankDeficientMatrix')
lso.RECT = true;  % linsolve opts, forces QR even when square (LU much worse)
co = linsolve(Q,d,lso);     % solve for proxy coeffs col vec
if o.verb, fprintf('cond(Q)=%.3g, norm(co)=%.3g\n',cond(Q),norm(co)); end

s8 = [];   % 3x3 near neighbors (excluding central) as source list...
for i=-1:1, for j=-1:1, if (i~=0 | j~=0)
      tr = i*e1+j*e2;        % translation in C
      s8 = [s8, bsxfun(@plus,s,[real(tr);imag(tr)])];
    end, end, end
if ich, ch8 = repmat(ch,[1 8]); else, ch8 = zeros(1,8*ns); end  % duplicate sources
if idip, dst8 = repmat(dst,[1 8]); dv8 = repmat(dv,[1 8]);
else, dst8 = []; dv8 = []; end
s8 = [s8 [real(p.x),imag(p.x)]'];     % append SLP proxy locs
ch8 = [ch8 (p.w.*co)'/cmcl]; dst8 = [dst8 zeros(1,M)]; dv8 = [dv8 zeros(2,M)];
ns8 = size(s8,2);

% *** case no self-evals for now (ipot,igr,ihe = 0)

% eval the original sources at self and/or targs...
O = lfmm2dpart(iprec,ns,s,ich,ch,idip,dst,dv,ipot,igr,ihe,nt,t,...
             ipott,igrt,ihet);
if ~o.noperi   % add in the periodizing part to answers
  % eval the 8 copies + proxies at targs...
  O8 = lfmm2dpart(iprec,ns8,s8,1,ch8,idip,dst8,dv8,0,0,0,nt,t,...
                  ipott,igrt,ihet);
  if ipott, O.pottarg = O.pottarg + O8.pottarg; end
  if igrt, O.gradtarg = O.gradtarg + O8.gradtarg; end
  if ihet, O.hesstarg = O.hesstarg + O8.hesstarg; end
end

function Q = Qmat(p,L,R,B,T,proxyrep) % matrix Q given proxy and colloc pts
% proxyrep is either LapSLP or LapDLP
[QL QLn] = proxyrep(L,p); [QR QRn] = proxyrep(R,p);
[QB QBn] = proxyrep(B,p); [QT QTn] = proxyrep(T,p);
Q = [QR-QL; QRn-QLn; QT-QB; QTn-QBn];

function f = discrepblk(f)    % sum 6 vector blocks with signs +++---
m = numel(f)/6; f = f(1:3*m)-f(3*m+1:6*m); f = f(1:m)+f(m+1:2*m)+f(2*m+1:end);


%%%%%%%
function test_lap2d2p    % tests and shows example usage
e1 = 1-.2i; e2 = 0.5+1i;   % lattice vecs as C numbers
iprec = 4;
rng(0);
ns = 10; s = e1*(rand(1,ns)-0.5)+e2*(rand(1,ns)-0.5);  % locs in UC
s = [real(s);imag(s)];  % 2-by-ns
ich = 0;ch = 0;   % no charge sources
idip=1; dz = randn(1,ns)+1i*randn(1,ns);   % random dipoles as complex #
dst = abs(dz); dv = [real(dz)./dst;imag(dz)./dst]; dv(isnan(dv)) = 0; % convert to strength, direction
ipot = 0; igr = 0; ihe = 0; % want no self-eval at srcs

ipott=1; igrt=1; ihet=0;        % periodizing accuracy test at 4 corners of UC
tt = -(e1+e2)/2 + [0 e1 e2 e1+e2]; tt = [real(tt);imag(tt)]; nt=size(tt,2);
o.verb = 0;
O = lap2d2p(iprec,ns,s,ich,ch,idip,dst,dv,ipot,igr,ihe,nt,tt,ipott,igrt,ihet,e1,e2,o);
u = real(O.pottarg); ue = max(u)-min(u);       % worst-case btw 4 corners
gu = real(O.gradtarg); gue = norm(max(gu,[],2)-min(gu,[],2));
fprintf('pointwise periodicity errors: potential %.3g, gradient %.3g\n',ue,gue)

ng = 100; x = 2*(1:ng)/ng-1; [xx yy] = meshgrid(x);  % fun plotting test
tt = [xx(:)';yy(:)']; nt = numel(xx);
ipott=1; igrt=0; ihet=0;
o.verb = 2;
[O U] = lap2d2p(iprec,ns,s,ich,ch,idip,dst,dv,ipot,igr,ihe,nt,tt,ipott,igrt,ihet,e1,e2,o);
u = real(O.pottarg);
imagesc(x,x,reshape(u,ng,ng)); caxis(3*[-1 1]+mean(u)); colorbar;
showunitcell(U); plot(s(1,:),s(2,:),'w.'); 
