function [crystalID1, crystalID2]=compute_crystalID_sub(transID1,transID2,axialID,transID,subset_id)
% compute crystal pair IDs for ordered subset using transaxial and axial IDs
%
% Author: Yansong Zhu @ UC Davis
% Last modified: 12/7/2023
num_sino_trans=length(transID1);
num_sino_ax=length(axialID);
num_crystal_trans=length(transID);
index_axial2=ceil(subset_id/(num_sino_trans*num_sino_ax));
tmp=subset_id-num_sino_trans*num_sino_ax*(index_axial2-1);
index_axial1=ceil(tmp/num_sino_trans);
index_trans=tmp-num_sino_trans*(index_axial1-1);

crystalID1=transID1(index_trans)+num_crystal_trans*axialID(index_axial1);
crystalID2=transID2(index_trans)+num_crystal_trans*axialID(index_axial2);
crystalID1=crystalID1(:);
crystalID2=crystalID2(:);
end