# PET-enabled dual-energy-CT
This code package provides a **real-scanner data implementation** for PET-enabled dual-energy CT (DECT). The code utilizes the projector modified from CASToR package, making this implementation work for general PET/CT scanners. This code package contains two different algorithms for the gamma-ray CT (gCT) image reconstruction from time-of-flight PET emission data: 
* Maximum likelihood of attenuation and activity (MLAA) 
* kernel MLAA.

A demo file is provided to show how to use this package to reconstruct gCT image from real scanner data. 

## Running the Demo:
To test the package, run Demo_testing_real_phantom.m.
Please download the data file from the following link:
https://drive.google.com/file/d/1LiUdmFmyugVr1URH0SjMHu71U7LSkIT6/view?usp=drive_link

## Use different PET scanners:
To adapt the code for other PET scanners, the following two modifications are required:
### Add scanner geometry file:
You will need to create a geometry file under /config/scanner/ describing your scanner’s geometry. You can use EXPLORER_oneUnit_24crystal.geom as a template. The geometry format follows the CASToR style. More details can be found in the [CASToR documentation](https://castor-project.org/sites/castor-project.org/files/2024-10/CASToR_general_documentation.pdf).
### Add system options in projector code:
Modify both proj_forw_CASTOR.m and proj_back_CASTOR.m to include your scanner’s input/output data settings. Use the case ‘EXPLORER_histogram’ as a reference when creating your own configuration.

## Note: 
1. Current gCT reconstruction code only supports sinogram format.
2. If you run the package on Linux system, the path for /config folder need to be added as environment varible so that the code could recognize the predefined PET system geometries.

If you use this package for your own research, please cite the corresponding research paper:\
[Y. Zhu, S. Li, Z. Xie, et al., “Feasibility of PET-enabled dual-energy CT imaging: First physical phantom and patient results”, European Journal of Nuclear Medicine and Molecular Imaging, 52, 1912-1923(2025).](https://link.springer.com/article/10.1007/s00259-024-06975-5)

## Other related work:
We also provide an open-source package for single-subject deep-learning gCT reconstruction, available here: https://github.com/SiqiLi1020/Single-subject-DL-reconstruction-with-optimization-transfer-for-PET-enabled-dual-energy-CT-imaging.


