# Supplementary Data Codes

**Three-dimensional imaging of individual carbon atoms**

Na Yeon Kim<sup>1*</sup>, Hanfeng Zhong<sup>1,2*</sup>, Jianhua Zhang<sup>3*</sup>, Colum M. O’Leary<sup>1*</sup>, Yuxuan Liao<sup>1</sup>, Ji Zou<sup>4</sup>, Haozhi Sha<sup>1</sup>, Minh Pham<sup>1</sup>, Weiyi Li<sup>1</sup>, Yakun Yuan<sup>1</sup>, Ji-Hoon Park<sup>5</sup>, Dennis Kim<sup>6</sup>, Huaidong Jiang<sup>3</sup>, Jing Kong<sup>5</sup>, Miaofang Chi<sup>7</sup>, Jianwei Miao<sup>1†</sup>,   

*<sup>1</sup>Department of Physics and Astronomy and California NanoSystems Institute, University of California, Los Angeles, CA, USA.*     
*<sup>2</sup>Department of Electrical and Computer Engineering, University of California, Los Angeles, CA, USA.*

*<sup>3</sup>Center for Transformative Science, ShanghaiTech University, Shanghai, China.*

*<sup>4</sup>Department of Physics, University of Basel, Basel, Switzerland.*

*<sup>5</sup>Department of Electrical Engineering and Computer Science, Massachusetts Institute of Technology, Cambridge, MA, USA.*

*<sup>6</sup>Department of Chemistry and Biochemistry, University of California, Los Angeles, CA, USA.*

*<sup>7</sup>Center for Nanophase Materials Sciences, Oak Ridge National Laboratory, Oak Ridge, TN, USA.*   

**These authors contributed equally to this work.*    

*†Correspondence author: J.M. (j.miao@ucla.edu).*     



## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Repositary Contents](#repositary-contents)
  
# Overview

# System Requirements

## Hardware Requirements

We recommend a computer with 32G DRAM, AMD Ryzen 7 6800H with Radeon Graphics 3.20 GHz CPU, and a NVIDIA GeForce RTX 3070 Ti GPU to run most data analysis source codes.

## Software Requirements

### OS Requirements

This package has been tested on the following Operating System:

Linux: Ubuntu 20.04.6 LTS (GNU/Linux 5.4.0-173-generic x86_64).   
Windows: Windows 11 Version 24H2 (OS Build 26100.3476).   
Mac OSX: We have not tested it on a Mac yet, but it should in principle work.   

### Matlab Version Requirements

This package has been tested with `Matlab` R2022b. All the codes have to run in their own folders. We recommend the use of `Matlab` version R2021a or higher to test the data and source codes.

# Repositary Contents

### 1. Experiment Data

Folder: [Measured_data](./1_Measured_data)

This folder contains 13 experimental projections after ptychography reconstruction and alignment as well as their corresponding angles.

### 2. The Extended Ptychographical Iterative Engine (EPIE) Package

Folder: [EPIE_package](./2_EPIE_package)

Run the ptychography reconstruction code `Main_EPIE.m` to obtain the phase retrieval projections of the twisted bilayer graphene sample.

### 3. The Probe REal Space Iterative REconstruction (RESIRE) Package

Folder: [RESIRE_package](./3_RESIRE_package)

Run the tomography reconstruction code `Main_AET_RESIRE_probe.m` to obtain the 3D reconstruction of the twisted bilayer graphene sample.

### 4. Reconstructed 3D Volume

Folder: [Final_reconstruction_volume](./4_Final_reconstruction_volume)

This folder contains the 3D volume of the twisted bilayer graphene reconstructed from `Main_AET_RESIRE_probe.m`.

### 5. Atom Tracing

Folder: [Tracing_atom_position](./5_Tracing_atom_position)

Run the code `Main_1_polynomial_tracing_iteration.m` to trace the initial atomic positions from the reconstructed 3D volume. After the manual checking of the 3D atomic positions, run the code `Main_2_match_atom_with_projection.m` to match the traced atom with the experimental projections. This is the pre-processing step for atimic position refinement.

### 6. Atomic Position Refinement

Folder: [Position_refinement](./6_Position_refinement)

Run the code `Main_atom_position_refinement.m` to refine the 3D atomic coordinates in the twisted bilayer graphene sample.

### 7. Experimental Atomic Model

Folder: [Final_coordinates](./7_Final_coordinates)

The final 3D atomic model and chemical element types (i.e. Carbon and Silicon) of the twisted bilayer graphene sample.

### 8. Match Experimental and Flat Atomic Models

Folder: [Match_atom](./8_Match_atom)

Run the code `Main_match_atom.m` to match the experimental atominc model with the flat atomic model.

### 9. Displacement Caculation

Folder: [Calculate_displacement](./9_Calculate_displacement)

Run the code `Main_calculate_displacement.m` to obtain the displacement between the experimental atomic model and the flat atomic model.

### 10. Meron and Skyrmion Caculation

Folder: [Calculate_Meron_Skyrmion](./10_Calculate_Meron_Skyrmion)

Run the 6 codes inside the folder `Main_calculate_antimeron_lower.m`,  `Main_calculate_antimeron_upper_1.m`,  `Main_calculate_antimeron_upper_2.m`, `Main_calculate_meron_lower.m`, `Main_calculate_meron_upper.m`, and  `Main_calculate_skyrmion_lower.m` to obatin the merons-like and skyrmion-like structure in the experimental atomic model as well as their corresponding skyrmion number.



