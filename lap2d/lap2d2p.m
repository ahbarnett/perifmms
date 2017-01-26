function [O U] = lap2d2p(iprec,ns,s,ich,ch,idip,dst,dv,ipot,igr,ihe,...
                         nt,t,ipott,igrt,ihet,e1,e2,o)
% draft of doubly periodic in 2D Laplace wrapper around CMCL FMMs.
%
% kernel is usual (1/2pi) log (1/r), which is -1/2pi times the odd CMCL
% normalization. Coordinates are as stacks of 2-cmpt real col vecs.
%
% unit cell centered at 0 with lattice vecs (e1,e2) each expressed as complex
% numbers. Last arg o is optional options struct.
% All other args, and output struct O, are as for lfmm2dpart.
%
% Barnett 1/26/17
if nargin==0, test_lap2d2p; return; end
if nargin<19, o=[]; end
if ~isfield(o,'verb'), o.verb=0; end
if ~isfield(o,'noperi'), o.noperi=0; end

badness = cond([real([e1 e2]);imag([e1 e2])]); % cond(lattice vector matrix)
% walls
U.e1 = e1; U.e2 = e2;
m = 25; [U L R B T] = doublywalls(U,m);  % note U is unit cell not CMCL U.
% proxy pts
proxyrep = @LapSLP;
M = ceil(80*badness);          % # proxy
Rp = 1.4; p.c = Rp * exp(1i*(0:M-1)'/M*2*pi);   % circle
p.x = e1*real(p.c) + e2*imag(p.c);     % ellipse
p = setupquad(p);         % adds s.w weights and s.nx normals to proxy
if o.verb, figure; showsegment({L R T B}); hold on; plot(p.x,'r.'); end
% periodizing solve Q block
warning('off','MATLAB:nearlySingularMatrix')  % backward-stable ill-cond is ok!
warning('off','MATLAB:rankDeficientMatrix')
lso.RECT = true;  % linsolve opts, forces QR even when square (LU much worse)
Q = Qmat(p,L,R,B,T,proxyrep);



% set up extended target list

% add near neighbors to source list


if o.noperi  % no periodization, for debug
  O=lfmm2dpart(iprec,ns,s,ich,ch,idip,dst,dv,ipot,igr,ihe,nt,t,...
               ipott,igrt,ihet);
end



function Q = Qmat(p,L,R,B,T,proxyrep) % matrix Q given proxy and colloc pts
% proxyrep is either LapSLP or LapDLP
[QL QLn] = proxyrep(L,p); [QR QRn] = proxyrep(R,p);
[QB QBn] = proxyrep(B,p); [QT QTn] = proxyrep(T,p);
Q = [QR-QL; QRn-QLn; QT-QB; QTn-QBn];

%%%%%%%
function test_lap2d2p
e1 = 1+.5i; e2 = -0.1+1i;
iprec = 4;
ns = 10; s = e1*(rand(1,ns)-0.5)+e2*(rand(1,ns)-0.5);
s = [real(s);imag(s)];  % 2-by-ns
ich = 0;ch = 0;   % no charge sources
idip=1; dz = randn(1,ns)+1i*randn(1,ns);   % random dipole as complex #
dst = abs(dz); dv = [real(dz)./dst;imag(dz)./dst]; dv(isnan(dv)) = 0; % convert to strength, direction
ipot = 0; igr = 0; ihe = 0; % want nothing at srcs
ng = 100; x = 2*(1:ng)/ng-1; [xx yy] = meshgrid(x);
tt = [xx(:)';yy(:)']; nt = numel(xx);
ipott = 1; igrt=0;ihet=0;
o.verb = 1; o.noperi = 1;
[O U] = lap2d2p(iprec,ns,s,ich,ch,idip,dst,dv,ipot,igr,ihe,nt,tt,ipott,igrt,ihet,e1,e2,o);
imagesc(x,x,reshape(real(O.pottarg),ng,ng)); caxis(ns*[-1 1]); colorbar;
showunitcell(U); plot(s(1,:),s(2,:),'w.'); 

% then check periodicity w/ 4 targs
