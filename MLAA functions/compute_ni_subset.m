function ni_sub=compute_ni_subset(crystal_eff_LUT,plane_eff,deadtime_LUT,trans_crystalID1,trans_crystalID2,sino_size_LR,sino_size_HR,subset_id)
% Compute multiplicative correction factors, including normalization
% factors and deadtime correction factors for subset LORs
%
% Inputs: 
%   crystal_eff_LUT: look-up table for crystal efficiency 
%   plane_eff:  look-up table for plane efficiency
%   deadtime_LUT: look-up table for deadtime correction factors
%   trans_crystalID1/trans_crystalID2: crystal ID pairs in transverse plane for each sinogram bin 
%   sino_size_LR: size for coarse sinogram (for scatter/random/deadtime)
%   sino_size_HR: size for fine sinogram (for emission data)
%   subset_id: IDs for subset sinogram bins
%
% Author: Yansong Zhu @ UC Davis
% Last modified: 11/28/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_crystal_trans,num_ax]=size(crystal_eff_LUT);
axialID=[1:num_ax]';
axialID=axialID-1;
num_sino_trans=length(trans_crystalID1);
index_axial2=ceil(subset_id/(num_sino_trans*num_ax));
tmp=subset_id-num_sino_trans*num_ax*(index_axial2-1);
index_axial1=ceil(tmp/num_sino_trans);
index_trans=tmp-num_sino_trans*(index_axial1-1);
%% for crystal eff
ID1=trans_crystalID1(index_trans(:))+num_crystal_trans*axialID(index_axial1(:))+1;
ID2=trans_crystalID2(index_trans(:))+num_crystal_trans*axialID(index_axial2(:))+1;
ni_sub_crystal=1./(crystal_eff_LUT(ID1).*crystal_eff_LUT(ID2)); ni_sub_crystal=ni_sub_crystal';
ni_sub_crystal=ni_sub_crystal(:);
clear ID1 ID2
%% for plane eff
ID1=axialID(index_axial1(:))+1;  
ID2=axialID(index_axial2(:));
ID=ID1+ID2*num_ax;
ni_sub_plaeff=1./plane_eff(ID);
clear ID1 ID2 ID
%% for deadtime
% compute axial ID for interpolation to match sinogram size with emission
% data
axial_id2=ceil(subset_id/(sino_size_HR(1)*sino_size_HR(2)));
tmp=subset_id-sino_size_HR(1)*sino_size_HR(2)*(axial_id2-1);
axial_id1=ceil(tmp/sino_size_HR(1));
trans_id=tmp-sino_size_HR(1)*(axial_id1-1);

axial1_grid=linspace(1,sino_size_HR(2),sino_size_LR(2));
axial2_grid=linspace(1,sino_size_HR(3),sino_size_LR(3));

if(length(axial1_grid)>1)
    axial1_diff=axial1_grid(2)-axial1_grid(1);
else
    axial1_diff=0;
end
if(length(axial2_grid)>1)
    axial2_diff=axial2_grid(2)-axial2_grid(1);
else
    axial2_diff=0;
end

axial1_left=discretize(axial_id1,axial1_grid);  axial1_right=axial1_left+1;  
axial2_left=discretize(axial_id2,axial2_grid);  axial2_right=axial2_left+1;  

LR_ID11=trans_id+sino_size_LR(1)*(axial1_left-1)+sino_size_LR(1)*sino_size_LR(2)*(axial2_left-1);
LR_ID12=trans_id+sino_size_LR(1)*(axial1_left-1)+sino_size_LR(1)*sino_size_LR(2)*(axial2_right-1);
LR_ID21=trans_id+sino_size_LR(1)*(axial1_right-1)+sino_size_LR(1)*sino_size_LR(2)*(axial2_left-1);
LR_ID22=trans_id+sino_size_LR(1)*(axial1_right-1)+sino_size_LR(1)*sino_size_LR(2)*(axial2_right-1);


x1=(axial1_grid(axial1_right)-axial_id1);
x2=(axial_id1-axial1_grid(axial1_left));
y1=(axial2_grid(axial2_right)-axial_id2);
y2=(axial_id2-axial2_grid(axial2_left));

clear axial1_right axial1_left axial2_right axial2_left axial1_grid axial2_grid

w11=x1.*y1; w12=x1.*y2; w21=x2.*y1; w22=x2.*y2;
DT_sub=repmat(w11,size(deadtime_LUT,1),1).*deadtime_LUT(:,LR_ID11)+...
       repmat(w12,size(deadtime_LUT,1),1).*deadtime_LUT(:,LR_ID12)+...
       repmat(w21,size(deadtime_LUT,1),1).*deadtime_LUT(:,LR_ID21)+...
       repmat(w22,size(deadtime_LUT,1),1).*deadtime_LUT(:,LR_ID22);
   
DT_sub=DT_sub/(axial1_diff*axial2_diff);    DT_sub=DT_sub';

ni_sub=ni_sub_crystal.*ni_sub_plaeff.*DT_sub;
end