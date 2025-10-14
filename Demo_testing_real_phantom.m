% This is a demo to test different image reconstruction algorithms for PET-eanbeld dual-energy CT on real phantom
% NOTE: Please do not distribute this package without permission from the authors

% For details of the reconstruction algorithms, please check the papers below:

% Standard MLAA:              A. Rezaei, M. Defrise, G. Bal, C. Michel, M. Conti, C. Watson, and J. Nuyts, "Simultaneous reconstruction of activity and attenuation in time-of-flight PET," IEEE Trans. Med.
%                             Imag., vol. 31, no. 12, pp. 2224-2233, Dec. 2012.

% Kernel MLAA (KAA):          G.B. Wang, "PET-enabled dual-energy CT: image reconstruction and a proof-of-concept computer simulation study," Phys. Med. Biol., vol. 65,
%                             pp. 245028, Nov. 2020.

% Programmer: Yansong Zhu, Siqi Li and Guobao Wang, UC DAVIS.
% Contact:    yszhu@health.ucdavis.edu; sqlli@health.ucdavis.edu; gbwang@health.ucdavis.edu
% Last modified:       10/14/2025


clc;
clear;
%% Load data
addpath('MLAA functions/');
load('Data/Sinogram_ScatterPhantom_24AxCrys_test.mat');
%% CT image with 80kvp as prior image
CT = CT_80kvp;
clear CT_image;
slice_per_unit = size(CT,3)/8;
CT_sub = CT(:,:,slice_per_unit*3+1:slice_per_unit*4);
CT_sub = CT_sub(:,:,26:42);
uinit = CT2LAC(CT_sub,'80','bilinear')/10;
%% quantification factor
QF = 1;
%% Setting
P = []; % system matrix for PET, required if mtype is set as 'matlab'
Popt.timefield = zeros(size(yi,2),1);
Popt.fovSize = [600,600,242.6/84*24];
Popt.num_event = size(yi,2);
% predefined PET scanner. To add a new scanner, please add scanner
% geometry file to /config/scanner, and other scanner parameters to 
% proj_forw_CASTOR.m/proj_back_CASTOR.m following the example of EXPLORER
% scanner
Popt.mtype = 'EXPLORER_histogram';
Popt.geo_file = 'EXPLORER_oneUnit_24crystal';
Popt.imgsiz = [150, 150, 17];
imgsiz = Popt.imgsiz;
Popt.disp = 1;
Popt.savestep = 1;
Popt.prjsiz = length(yi(:))/27;

%% for HR data
Popt.transID1 = crystalID1_trans_HR'; % transaxial crystal ID pairs for each LOR in sinogram
Popt.transID2 = crystalID2_trans_HR';
Popt.axID = [0:23]';
Popt.transID_table = [1:840]';
% the following three lines are used to compute normalization factors
Popt.crys_eff_sino_LUT = cryseff;
Popt.DT_sino = DT_sino_reciprocal;
Popt.planeff = planeff;
%
Popt.nbTOFBins = 27; % number of TOF bins in sinogram
% low-resolution sinogram size used for scatter/random/deadtime data(numbin*numproj,axialID1,axialID2)
Popt.sino_size_LR = [533*420,2,2]; 
% high-resolution sinogram size for emission data yi
Popt.sino_size_HR = [533*420,24,24];

%% random/scatter events
% in this example, random/scatter events are stored in a low-resolution sinogram
% to save storage/memory usage. Interpolation will be performed during
% reconstruction to match the size of low-resolution/high-resolution
% sinogram
ri = rand_event + scat_event;
clear rand_event scat_event

%% MLAA setting
Gopt_attn = Popt;
Gopt_attn.imgsiz = [150, 150, 17];
Gopt_attn.mtype = 'EXPLORER_histogram';
Gopt_attn.geo_file = 'EXPLORER_oneUnit_24crystal';
G_attn = []; % system matrix for CT, required if mtype is set as 'matlab'

%% Set TOF flag for Popt/Gopt
Popt.TOF_flag = 1;
Gopt_attn.TOF_flag = 0;

%% Algorithm setting
xinit = 30*ones(imgsiz);
Gopt_attn.mask = uinit > 0;
Popt.mask = xinit > 0;
% Total iteration = maixt * number of subset
maxit = 20;
num_subset = 20;

%% build kernel (uncomment this part if using kernel MLAA method)
% fprintf('---Start building kernel---\n')
% imgsiz_CT = size(CT_sub);
% R = buildNbhd(imgsiz_CT, 'clique', 1); % Extract features using a 3x3 patch 
% I = [[1:prod(imgsiz_CT)]' R.N]; 
% F = CT_sub(I);
% F = F * diag(1./std(F,1)); % normalization 
% % building kernel matrix K 
% sigma = 1;
% [N, W] = buildKernel(imgsiz, 'knn', 50, F, 'radial', sigma);
% K = buildSparseK(N, W);
%% Choose reconstrution alogrithm
% MLAA, KAA
rectype = 'MLAA';

%% Synergistic reconstruction

switch rectype

    case 'MLAA'
        disp('------ Standard MLAA reconstruction...')
        tic;
        [u, x, out] = psct_kmlaa_CASTOR_OS(yi, G_attn, Gopt_attn, uinit, P, Popt, xinit, ri, maxit, num_subset);
        toc;
        out.uest_withoutK=out.uest;
        u=double(u);
        toc;
        out_store=out;
        currentMaxIt = maxit * num_subset;
        save('testing_result/MLAA_real_phantom.mat','num_subset','xinit','uinit','out_store','currentMaxIt','-v7.3')

    case 'KAA'
            disp('------ Kernel MLAA reconstruction...')
            tic;
            [u, x, out] = psct_kmlaa_CASTOR_OS(yi, G_attn, Gopt_attn, uinit, P, Popt, xinit, ri, maxit, num_subset, K);
            toc;
            out.uest_withoutK=out.uest;
            out.uest =K*double(out.uest);
            u=double(u);
            out_store=out;
            currentMaxIt = maxit * num_subset;
            save('testing_result/KAA_real_phantom.mat','num_subset','xinit','uinit','out_store','currentMaxIt','-v7.3')
end
