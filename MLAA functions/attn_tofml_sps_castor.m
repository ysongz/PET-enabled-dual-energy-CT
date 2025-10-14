function [x, li, L] = attn_tofml_sps_castor(yt, bt, G, Gopt, x, rt, maxit, li, aa)
% This is a sub program for attenuation map reconstruction from PET emssion
% data
%
% gbwang@ucdavis.edu 05-01-2018
%

%% check inputs

numprj = prod(Gopt.prjsiz);
numfrm = size(yt,2);
numbin = size(yt,1)/numprj;

if nargin<2 | isempty(bt)
    bt = ones(numprj,numfrm);
end
if nargin<6 | isempty(rt)
    rt = ones(numprj*numbin,numfrm);
    yt = yt + rt;
end
if nargin<8
    li = [];
end
if nargin<9 
    aa = [];
end
if isempty(aa)
    aa = proj_forw_CASTOR(G, Gopt, ones(size(x)));
end
%% convert data precision
yt=single(yt);  bt=single(bt);  G=single(G);
x=single(x);    rt=single(rt);  li=single(li);
aa=single(aa);

rt=rt+eps;
%% iterate
for it = 1:maxit
    
    % foward projection
    if isempty(li)
        li = proj_forw_CASTOR(G, Gopt, x);
    end
    
    
    % gradient
    lt = li'; lt=repmat(lt,[numbin 1]); % repeat for TOF bins
    lt = repmat(lt(:),[1 numfrm]); % repeat for number of frames
    
    yb = bt.*exp(-lt) + rt;
    yr = (1-yt./yb).*(yb-rt);
    yr = reshape(yr, [numbin numprj numfrm]);
    hi = sum(sum(yr,3),1); hi=hi';
    gx = proj_back_CASTOR(G, Gopt, hi);
    clear hi yr
    % optimal curvature and the optimization transfer weight
    nt = sum(trl_curvature_castor(yt, bt, rt, lt, 'oc'),2);
    nt = reshape(nt,[numbin numprj]);   nt=nt';
    wx = proj_back_CASTOR(G, Gopt, sum(nt,2).*aa);% nt has been summed over nbin and nframe
    
    % objective function
    L(it) = sum(yt(:).*log(yb(:)) - yb(:));
    
    % update image
    x(Gopt.mask) = x(Gopt.mask) + gx(Gopt.mask) ./ wx(Gopt.mask);
    x(Gopt.mask(:)&wx(:)==0) = 0;
    x = max(0,x);
    
    % update projection
    if it<maxit
        li = proj_forw_CASTOR(G, Gopt, x, numfrm);
    end
    
end
