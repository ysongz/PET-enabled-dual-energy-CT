function ri_sub=compute_ri_subset(ri,subset_id,sino_size_LR,sino_size_HR)
% Compute additive correction factors, including scatter and random factors for subset LORs
% Interpolation is performed to match sinogram size of ri and emission data
% 
% Inputs:
%   ri: sinogram of additive factors
%   subset_id: IDs for subset sinogram bins
%   sino_size_LR: size for coarse sinogram (for scatter/random/deadtime)
%   sino_size_HR: size for fine sinogram (for emission data)
%
% Author: Yansong Zhu @ UC Davis
% Last modified: 11/28/2023
%% compute axial ID for interpolation
axial_id2=ceil(subset_id/(sino_size_HR(1)*sino_size_HR(2)));
tmp=subset_id-sino_size_HR(1)*sino_size_HR(2)*(axial_id2-1);
axial_id1=ceil(tmp/sino_size_HR(1));
trans_id=tmp-sino_size_HR(1)*(axial_id1-1);

% trans_grid=linspace(1,sino_size_HR(1),sino_size_LR(1));
axial1_grid=linspace(1,sino_size_HR(2),sino_size_LR(2));
axial2_grid=linspace(1,sino_size_HR(3),sino_size_LR(3));
% trans_diff=trans_grid(2)-trans_grid(1);
axial1_diff=axial1_grid(2)-axial1_grid(1);
axial2_diff=axial2_grid(2)-axial2_grid(1);

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

w11=x1.*y1; w12=x1.*y2; w21=x2.*y1; w22=x2.*y2;
ri_sub=repmat(w11,size(ri,1),1).*ri(:,LR_ID11)+...
       repmat(w12,size(ri,1),1).*ri(:,LR_ID12)+...
       repmat(w21,size(ri,1),1).*ri(:,LR_ID21)+...
       repmat(w22,size(ri,1),1).*ri(:,LR_ID22);
   
ri_sub=ri_sub/(axial1_diff*axial2_diff);
%% normalize for different sinogram size
%sf=sino_size_LR(2)*sino_size_LR(3)/sino_size_HR(2)/sino_size_HR(3);
%ri_sub=ri_sub*sf;
end
