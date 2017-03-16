function [O U] = lfmm2d2ppart(iprec,ns,s,ich,ch,idip,dst,dv,ipot,igr,ihe,...
                         nt,t,ipott,igrt,ihet,v1,v2,o)
% lfmm2d2ppart  Doubly-periodic FMM for 2D laplace kernel, general unit cell
%
% [O U] = lfmm2d2ppart(iprec,ns,s,ich,ch,idip,dst,dv,ipot,igrad,ihess,...
%                         nt,t,ipott,igradt,ihesst,v1,v2,opts)
% Inputs:
%  First 16 inputs (iprec,....,ihesst): same as inputs to CMCL lfmm2dpart
%  v1, v2 : 2x1 real-valued vectors giving lattice vectors
%           Note the orientation sense: angle(v2,v1)<pi, ie cross(v1,v2)>0
%  opts   : optional structure with optional fields:
%           opts.verb = 0 (silent), 1 (text output), 2 (text & diagnostic figs)
%           opts.noperi : if true, don't periodize
% Outputs:
%  O      : same as outputs of CMCL lfmm2dpart
%  U      : unit cell struct with fields giving wall points, etc
%
% This evaluates the 2D Laplace potential (and gradient, hessian, if requested)
% due to the given set of sources (and dipoles, if requested) doubly-periodized
% by the lattice with vectors v1 and v2. It uses the CMCL FMMLIB2D code, and
% has a nearly identical interface (however, see below re the different
% normalization). This means for N sources and N targets, the
% run time is O(N). It can handle evaluation at the sources, in which case
% sum is over terms i.neq.j, plus additional targets (for which the sum has
% no such restriction). The resulting potential is only defined up to a
% constant.
%
% Assumptions:
%  * charges sum to zero, for compatibility w/ periodizing.
%  * all sources and targets lie within or close to the unit cell
%  * unit cell is centered on origin (seems reasonable)
%
% The kernel is the usual (1/2pi) log (1/r). Note this is -1/2pi times the CMCL
% normalization. Coordinates are as stacks of 2-cmpt real col vecs.
% Source and dipole strengths are real valued.
%
% See also: lfmm2dpart from CMCL FMMLIB2D

% todo:
% 1) think about reading off and enforcing nonzero potential drops R-L and
%  T-B from the requested sources...
% 2) more tests including wacky large/small/skew unit cells.
%
% Alex Barnett
%   initial version afternoon of 1/26/17
%   additional non-src targets 1/28/17
%   proxy ellipse and 3x3 -> proxy box and 2x2; around twice speed. 3/5/17
%   real not complex unit cell vec interface, name change, rel errs. 3/8/17
if nargin==0, test_lfmm2d2ppart; return; end
if nargin<19, o=[]; end
if ~isfield(o,'verb'), o.verb=0; end           % how much diagnostic output
if ~isfield(o,'noperi'), o.noperi=0; end
e1 = v1(1)+1i*v1(2); e2 = v2(1)+1i*v2(2);      % unit cell vecs as C#s.
AU = [v1 v2];          % lattice vector matrix
badness = max(1.3,cond(AU));    % cond(lattice vector matrix)
if badness>10, warning('unit cell close to singular: will be slow & bad!'); end
cmcl = -2*pi;          % factor by which CMCL normalization is off
ch = ch(:)'*(1/cmcl); dst = dst(:)'*(1/cmcl);  % correct & ensure row vecs
digitlist = [0 1 2 3 6 9 12 15]; digits = digitlist(iprec+3);  % log10(tol)

% collocation walls
tt=tic;             % internal timings
U.e1 = e1; U.e2 = e2; U.v1 = v1; U.v2 = v2;
m = ceil(5+1.3*digits*badness);            % # colloc pts per wall
[U L R B T] = doublywalls(U,m);  % note U is unit cell not CMCL U.
n1 = T.nx(1); n2 = R.nx(1);           % the two wall normals
% proxy pts, uniform on (fixed) 2x scaled unit cell
M = ceil(5+1.6*digits*badness);          % # proxy per wall
z = 2*(0.5:M-0.5)'/M-1;            % col vec of nodes in [-1,1)
p.x = [-e2+z*e1; e1+z*e2; e2+z*e1; -e1+z*e2];  % box
p.nx = kron([-n1;n2;n1;-n2],ones(M,1));  % outward normals (only needs for DLP)
p.w = 0*p.x+2/M;              % col vec of weights of order the spacing
if o.verb>1, figure; showsegment({L R T B}); hold on; plot(p.x,'r.'); end
% periodizing Q block
Q = Qmat(p,L,R,B,T,@LapSLP);    % use proxy monopoles
if o.verb, fprintf('Q fill\t\t%.3f s;\tbadness %.3g, Nprox=%d, m=%d\n',toc(tt),badness,4*M,m), end

% compute discrep: one src to 8 targ wall copies (ns -> around 200 targs)
tt=tic;
tw = [R.x-e2;R.x;L.x-e1-e2;L.x-e1; T.x-e1;T.x;B.x-e1-e2;B.x-e2].'; %row vec in C
tw = [real(tw);imag(tw)];   % 2-by-n format
s0 = inv(AU)*s;             % srcs mapped back to [-1/2,1/2]^2
if sum(abs(s0(:))>0.6), warning('source points further than 0.1 outside unit cell: probably inaccurate!'); end
sd = s - v1*(s0(1,:)>0) - v2*(s0(2,:)>0);  % shift down 1 cell the >0 coords
W = l2dpartdirect(ns,sd,ich,ch,idip,dst,dv,0,0,0,8*m,tw,1,1,0); % beats fmm
f = discrepblk(W.pottarg(1:4*m));       % 1st discrep block
fn = discrepblk([real(n2),imag(n2)] * W.gradtarg(:,1:4*m));  % wall n-derivs
g = discrepblk(W.pottarg(4*m+1:end));
gn = discrepblk([real(n1),imag(n1)] * W.gradtarg(:,4*m+1:end));
d = -real([f,fn,g,gn])';     % cancel out discrep col vec
if o.verb, fprintf('discrep\t\t%.3f s\n',toc(tt)), end

tt=tic;
warning('off','MATLAB:nearlySingularMatrix')  % backward-stable ill-cond is ok!
warning('off','MATLAB:rankDeficientMatrix')
lso.RECT = true;  % linsolve opts, forces QR even when square (LU much worse)
co = linsolve(Q,d,lso);     % solve for proxy coeffs col vec
if o.verb, fprintf('linsolve\t%.3f s;\tnorm(co)=%.3g\n',toc(tt),norm(co)); end

% eval the original sources at self and/or targs...
tt=tic;
O = lfmm2dpart(iprec,ns,s,ich,ch,idip,dst,dv,ipot,igr,ihe,nt,t,...
             ipott,igrt,ihet);
if o.verb, fprintf('self FMM\t%.3f s\n',toc(tt)), end

if ~o.noperi   % add in the periodizing part to answers
  s3 = [];   % 2x2 near neighbors (excluding centered 1x1 block) as source list...
  t1s = v1*(2*(s0(1,:)>0)-1);   % sgn(x)*v1 for each pt (x,y) back in [-1/2,1/2]^2
  t2s = v2*(2*(s0(2,:)>0)-1);   % sgn(y)*v2
  for i=-1:0, for j=-1:0, if (i~=0 | j~=0)   % 3 translations, ie L shape
        s3 = [s3, s + i*t1s + j*t2s];
      end, end, end
  if ich, ch3 = repmat(ch,[1 3]); else, ch3 = zeros(1,3*ns); end  % duplicate srcs
  if idip, dst3 = repmat(dst,[1 3]); dv3 = repmat(dv,[1 3]);
  else, dst3 = []; dv3 = []; end
  s3 = [s3 [real(p.x),imag(p.x)]'];     % append SLP proxy locs
  ch3 = [ch3 (p.w.*co)'/cmcl]; dst3 = [dst3 zeros(1,4*M)]; dv3 = [dv3 zeros(2,4*M)];
  ns3 = size(s3,2);
      
  srceval = ipot || igr || ihe; % if so, append to targs list for periodizing
  if srceval, t = [s t]; nst = ns+nt; ti=ns+(1:nt);   % ti = targ inds in t
    ipott3 = ipot | ipott; igrt3 = igr | igrt; ihet3 = ihe | ihet;
  else
    nst = nt; ti = 1:nt; ipott3=ipott; igrt3=igrt; ihet3=ihet;  % just targs
  end
  % eval the 3 copies + proxies at all targs (of course no self-eval here)...
  tt=tic;
  O3 = lfmm2dpart(iprec,ns3,s3,1,ch3,idip,dst3,dv3,0,0,0,nst,t,...
                  ipott3,igrt3,ihet3);
  if o.verb, fprintf('3N nei+pxy FMM\t%.3f s\n',toc(tt)), end
  if ipot, O.pot = O.pot + O3.pottarg(1:ns); end
  if igr, O.grad = O.grad + O3.gradtarg(:,1:ns); end
  if ihe, O.hess = O.hess + O3.hesstarg(:,1:ns); end
  if ipott, O.pottarg = O.pottarg + O3.pottarg(ti); end
  if igrt, O.gradtarg = O.gradtarg + O3.gradtarg(:,ti); end
  if ihet, O.hesstarg = O.hesstarg + O3.hesstarg(:,ti); end
end

function Q = Qmat(p,L,R,B,T,proxyrep) % matrix Q given proxy and colloc pts
% proxyrep is either LapSLP or LapDLP
[QL QLn] = proxyrep(L,p); [QR QRn] = proxyrep(R,p);
[QB QBn] = proxyrep(B,p); [QT QTn] = proxyrep(T,p);
Q = [QR-QL; QRn-QLn; QT-QB; QTn-QBn];

function f = discrepblk(f)    % sum 4 vector blocks with signs ++--
m = numel(f)/4; f = f(1:2*m)-f(2*m+1:4*m); f = f(1:m)+f(m+1:end);


%%%%%%%
function test_lfmm2d2ppart     % tests and shows example usage
v1 = [1; -.2]; v2 = [.5; 1];   % lattice vectors (must be column vectors)
%v1 = [1; 1]; v2 = [.5; 1];     % badness 6
%v1 = [1;0]; v2 = [0;.1];       % aspectratio 10
%v1 = [.1;0]; v2 = [0;1];       % "
rng(0);    % fix random seed
ns = 1e4;   % number of sources ...try 10 or 1e4 or more
s = v1*(rand(1,ns)-0.5)+v2*(rand(1,ns)-0.5);  % 2*ns locs in UC
ich = 0; ch = 0;   % no charge sources
idip=1; dz = randn(1,ns)+1i*randn(1,ns);   % random dipoles as complex #
dst = abs(dz); dv = [real(dz)./dst;imag(dz)./dst]; dv(isnan(dv)) = 0; % convert to strength, direction
ipot = 1; igr = 1; ihe = 0; % what is wanted for self-interactions (i.neq.j)

ipott=1; igrt=1; ihet=0;        % periodizing accuracy test at 4 corners of UC
e1 = v1(1)+1i*v1(2); e2 = v2(1)+1i*v2(2);      % unit cell vecs as C#s
tt = -(e1+e2)/2 + [0 e1 e2 e1+e2]; tt = [real(tt);imag(tt)]; nt=size(tt,2);
o.verb = 0;
for iprec=1:5   % check different requested accuracies
  O = lfmm2d2ppart(iprec,ns,s,ich,ch,idip,dst,dv,ipot,igr,ihe,nt,tt,ipott,igrt,ihet,v1,v2,o);
  % extract potential and gradient, get their errors...
  u = real(O.pottarg); ue = max(u)-min(u);    % worst-case btw 4 corners
  gu = real(O.gradtarg); gue = norm(max(gu,[],2)-min(gu,[],2));  % "
  uer = ue/max(abs(u)); guer = gue/norm(max(gu,[],2));    % relative errors
  fprintf('iprec=%d: ptwise periodicity relative errors: pot %.3g, grad %.3g\n',iprec,uer,guer)
end
%u = real(O.pot(2:end) - O.pot(1))  % check self-eval converges varying M, etc
%real(O.grad)   % check conv

iprec = 4;
ng = 100; x = 2*(1:ng)/ng-1; [xx yy] = meshgrid(x);  % fun plotting test
tt = [xx(:)';yy(:)']; nt = numel(xx);
ipott=1; igrt=0; ihet=0;
o.verb = 1; tic
[O U] = lfmm2d2ppart(iprec,ns,s,ich,ch,idip,dst,dv,ipot,igr,ihe,nt,tt,ipott,igrt,ihet,v1,v2,o);
fprintf('%d periodized src to %g targs in %.3g s\n',ns,nt,toc)
u = reshape(real(O.pottarg),ng,ng);   % extract potential
u0 = u(round(ng/2),round(ng/2));      % central value for colors
imagesc(x,x,u); caxis(1.5*sqrt(ns)*[-1 1]+u0); colorbar; hold on;
showunitcell(U); if ns<1e3, plot(s(1,:),s(2,:),'w.'); end
