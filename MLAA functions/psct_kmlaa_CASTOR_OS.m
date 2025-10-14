function [u, x, out] = psct_kmlaa_CASTOR_OS(yi, A, Aopt, u, G, Gopt, x, ri, maxit, num_subset,K_CT,id_subset)
%
% CT-kernelized MLAA algorithm for synergistic reconstructon of time-of-flight 
% PET activity and attenuation images
%
% gbwang@ucdavis.edu Aug 02, 2019
%
% Last updated: March 12, 2021
%
% modified by Yansong, 12/1/2021
% ordered subset version for kMLAA/MLAA
%% check inputs for reconstruction

% total frame number
numprj = prod(Aopt.prjsiz);
numbin = size(yi,1);
% CT operator
Aopt = setGopt_castor(ones(numprj,1), A, Aopt);
numpix_CT = prod(Aopt.imgsiz);
if nargin<5 | isempty(u)
    u = zeros(numpix_CT,1);
end

% PET operator
numpix = prod(Gopt.imgsiz);
Gopt = setGopt_castor([], G, Gopt);
if ~isfield(Gopt,'emisFlag')
    Gopt.emisFlag = 1;
end
if ~isfield(Aopt,'attnFlag')
    Aopt.attnFlag = 1;
end
if nargin<8 | isempty(x)
    x = ones(numpix,1);
end
if Gopt.emisFlag
    x = max(mean(x(:))*1e-9,x(:)); 
    x(~Gopt.mask,:) = 0;
end
if nargin<8 | isempty(ri)
    yeps = mean(yi(:))*1e-9;
    ri = ones(size(yi))*yeps;
    yi = yi + ri;
end

% iteration number
if nargin<9 | isempty(maxit)
    maxit = 10;
end

% kernel for gamma-ray CT reconstruction
if nargin<11 | isempty(K_CT)
    K_CT = speye(numpix_CT);
end
Aopt.kernel = K_CT;
%% convert data precision  
x=single(x);    
u=single(u);
% iteration number
if isfield(Aopt,'MaxIt')
    it_attn = Aopt.MaxIt;
else
    it_attn = 1;
end
if isfield(Gopt,'MaxIt')
    it_emis = Gopt.MaxIt;
else
    it_emis = 1;
end



%% For ordered subset
% reshape yi, ni, ri
num_TOFbin=length(yi)/Gopt.num_event;

if nargin<12 | isempty(id_subset) % if subset sequence is not provided as input, generate random subset
% random index for ordered subset
    id=randperm(Gopt.num_event);
    id_sub=cell(num_subset,1);
    event_persubset=round(Gopt.num_event/num_subset);
    for i=1:num_subset-1
        id_sub{i}=id((i-1)*event_persubset+1:i*event_persubset);
    end
    id_sub{num_subset}=id((num_subset-1)*event_persubset+1:end);
    clear id
else
    id_sub=cell(num_subset,1);
    if(length(id_subset)~=num_subset)
        disp('number of subset does not match subset sequence');
    else
        for i=1:num_subset
            id_sub{i}=id_subset{i};
        end
    end
end

if ~Aopt.attnFlag
    gg=cell(num_subset,1);
    for i=1:num_subset
       Atmp=A;
       Aopt_tmp=subset4opt(Aopt,id_sub{i});
       Gtmp=G;
       Gopt_tmp=subset4opt(Gopt,id_sub{i});
       li_attn = proj_forw_CASTOR(Atmp, Aopt_tmp, u);
       ai = exp(-li_attn);
       ai=ai'; ai=repmat(ai,[numbin 1]);
       ni_sub=compute_ni_subset(Gopt_tmp.crys_eff_sino_LUT,Gopt_tmp.planeff,Gopt_tmp.DT_sino,...
                                 Gopt_tmp.transID1,Gopt_tmp.transID2,...
                                 Gopt_tmp.sino_size_LR,Gopt_tmp.sino_size_HR,id_sub{i});
       ni_sub=ni_sub'; ni_sub=repmat(ni_sub,[numbin 1]);
       mi = ni_sub(:).*ai(:);
       gg{i} = proj_back_CASTOR(G, Gopt, mi(:));
    end
else
    aa=cell(num_subset,1);
    for i=1:num_subset
        Atmp=A;
        Aopt_tmp=subset4opt(Aopt,id_sub{i});
        aa{i} = proj_forw_CASTOR(Atmp, Aopt_tmp, ones(numpix_CT,1)); % aa = A*1
    end
end
if ~Gopt.emisFlag
    li_emis_store=cell(num_subset,1);
    for i=1:num_subset
        Gopt_tmp=subset4opt(Gopt,id_sub{i});
        li_emis_store{i} = proj_forw_CASTOR(G, Gopt_tmp, x);
    end
end

% output
if nargin>1
    out = []; L = [];
end
out.xest = zeros(length(x(:)), ceil(maxit/Gopt.savestep)+1);
out.uest = zeros(length(u(:)), ceil(maxit/Gopt.savestep)+1);
t1 = tic;
%% iterative loop
for it = 1:maxit     
    disp(it);
    % save data
    if nargout>1 & ( it==1 | rem(it,Gopt.savestep)==0 )
        itt = min(it, floor(it/Gopt.savestep) + 1);
        if Gopt.disp==1
            disp(sprintf('iteration %d',it));
        end
        out.step(itt)   = it-1;
        out.time(:,itt) = toc(t1);        
        out.xest(:,itt) = x(:);
        out.uest(:,itt) = u(:);
    end
	for it_sub=1:num_subset
		disp(sprintf('subset iteration %d',it_sub))
        Gtmp=G;
        Gopt_tmp=subset4opt(Gopt,id_sub{it_sub});
        Atmp=A;
        Aopt_tmp=subset4opt(Aopt,id_sub{it_sub});
        ni_sub=compute_ni_subset(Gopt.crys_eff_sino_LUT,Gopt.planeff,Gopt.DT_sino,...
                              Gopt.transID1,Gopt.transID2,...
                              Gopt.sino_size_LR,Gopt.sino_size_HR,id_sub{it_sub});
        ni_sub=repmat(ni_sub',Gopt.nbTOFBins,1);
        ri_sub=compute_ri_subset(ri,id_sub{it_sub},Gopt.sino_size_LR,Gopt.sino_size_HR);
        ri_sub=single(ri_sub); ri_sub=ri_sub.*ni_sub;% needed if ri_sub is normalization-corrected

        yi_sub=yi(:,id_sub{it_sub});    yi_sub=yi_sub(:); yi_sub=full(yi_sub);
        yi_sub=single(yi_sub);
        if(it_sub<num_subset)
            Aopt_tmpnext=subset4opt(Aopt,id_sub{it_sub+1});
            Gopt_tmpnext=subset4opt(Gopt,id_sub{it_sub+1});
        else
            Aopt_tmpnext=subset4opt(Aopt,id_sub{1});
            Gopt_tmpnext=subset4opt(Gopt,id_sub{1});
        end
        if Gopt.emisFlag && it==1 && it_sub==1
            li_emis = proj_forw_CASTOR(Gtmp, Gopt_tmp, x);
        end
        if Aopt.attnFlag && it==1 && it_sub==1
            li_attn = proj_forw_CASTOR(Atmp, Aopt_tmp, u);
        end
        if ~Gopt.emisFlag
            li_emis=li_emis_store{it_sub};
        end
        % likelihood function
        ai=exp(-li_attn);
        ai=ai';	ai=repmat(ai,[numbin,1]);
        yb = ni_sub(:).*ai(:).*li_emis(:) + ri_sub(:);
        L(it) = sum(yi_sub(:).*log(yb(:)+eps)-(yb(:)+eps));
    
        % attenuation estimate
        if Aopt.attnFlag 
            tic
            mi = ni_sub(:).*li_emis(:);
            [u, V, Lattn] = attn_tofml_sps_castor(yi_sub, mi(:), Atmp, Aopt_tmp, u(:), ri_sub(:), it_attn, li_attn(:), aa{it_sub});
            li_attn = proj_forw_CASTOR(Atmp, Aopt_tmp, u);
            li_attn_next=proj_forw_CASTOR(Atmp, Aopt_tmpnext, u);
            ai=exp(-li_attn);
            ai=ai'; ai=repmat(ai,[numbin 1]); 
            toc
            L_attn(it)=Lattn;
        end
    
    % emission estimate
        if Gopt.emisFlag
            tic
            if Aopt.attnFlag 
                mi = ni_sub(:).*ai(:);
                gg{it_sub} = proj_back_CASTOR(Gtmp, Gopt_tmp, mi);
            end
            [x, U] = emis_ml_em_castor(yi_sub(:), mi(:), Gtmp, Gopt_tmp, x(:), ri_sub(:), it_emis, li_emis(:), gg{it_sub});   
            li_emis = proj_forw_CASTOR(Gtmp, Gopt_tmpnext, x);
            toc
        end
        li_attn=li_attn_next;
    end
end
% save final iteration
itt=ceil(maxit/Gopt.savestep)+1;
if Gopt.disp==1
   disp(sprintf('save final iteration'));
end
out.step(itt)   = it;  
out.time(:,itt) = toc(t1);
out.xest(:,itt) = x(:);
out.uest(:,itt) = u(:);
% output
out.cost = L;

proj_clear(Gopt);
end

function Gopt_sub=subset4opt(Gopt,id_sub)
%% set system information for subset
Gopt_sub=Gopt;
[tmpID1,tmpID2]=compute_crystalID_sub(Gopt.transID1,Gopt.transID2,Gopt.axID,Gopt.transID_table,id_sub);
Gopt_sub.crystalID1=tmpID1;
Gopt_sub.crystalID2=tmpID2;
clear tmpID1 tmpID2
Gopt_sub.timefield=Gopt.timefield(id_sub);
Gopt_sub.num_event=length(id_sub);
Gopt_sub.prjsiz=length(id_sub);
end
