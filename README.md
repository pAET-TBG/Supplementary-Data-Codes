# Supplementary Data Codes

# **Three-dimensional determination of individual carbon atoms**

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

This repository provides all experimental and simulated data, along with source code for ptychographic and tomographic reconstruction, atom tracing, position refinement, and data analysis associated with our study. We report the first experimental determination of the 3D coordinates of individual carbon atoms with picometer-scale precision using ptychographic atomic electron tomography (pAET). Using twisted bilayer graphene as a model system, our results reveal chiral lattice distortions and interlayer interactions that were previously inaccessible to conventional microscopy. We further demonstrate the generalizability of pAET to disordered and low-Z materials through simulations on amorphous carbon and a nitrogen-vacancy center in diamond. All resources are freely available for download and unrestricted use by the research community. We hope this repository supports further advances in high-precision 3D atomic imaging.

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

### a. Experiment Data

Folder: [Measured_data](./a_Measured_data)

This folder contains 13 experimental projections after ptychography reconstruction and alignment as well as their corresponding angles.

### b. The Extended Ptychographical Iterative Engine (EPIE) Package

Folder: [EPIE_package](./b_EPIE_package)

Run the ptychography reconstruction code `Main_EPIE.m` to obtain the phase retrieval projections of the twisted bilayer graphene sample. Run the codes inside Main_prepare_tomography to obtain the projections and angles for tomography

### c. The Probe REal Space Iterative REconstruction (RESIRE) Package

Folder: [RESIRE_package](./c_RESIRE_package)

Run the tomography reconstruction code `Main_AET_RESIRE_probe.m` to obtain the 3D reconstruction of the twisted bilayer graphene sample.

### d. Reconstructed 3D Volume

Folder: [Final_reconstruction_volume](./d_Final_reconstruction_volume)

This folder contains the 3D volume of the twisted bilayer graphene reconstructed from `Main_AET_RESIRE_probe.m`.

### e. Atom Tracing

Folder: [Tracing_atom_position](./e_Tracing_atom_position)

Run the code `Main_1_polynomial_tracing_iteration.m` to trace the initial atomic positions from the reconstructed 3D volume. After the manual checking of the 3D atomic positions, run the code `Main_2_match_atom_with_projection.m` to match the traced atom with the experimental projections. This is the pre-processing step for atimic position refinement.

### f. Atomic Position Refinement

Folder: [Position_refinement](./f_Position_refinement)

Run the code `Main_atom_position_refinement.m` to refine the 3D atomic coordinates in the twisted bilayer graphene sample.

### g. Experimental Atomic Model

Folder: [Final_coordinates](./g_Final_coordinates)

The final 3D atomic model and chemical element types (i.e. Carbon and Silicon) of the twisted bilayer graphene sample.

### h. Match Experimental and Flat Atomic Models

Folder: [Match_atom](./h_Match_atom)

Run the code `Main_match_atom.m` to match the experimental atominc model with the flat atomic model.

### i. Displacement Caculation

Folder: [Calculate_displacement](./i_Calculate_displacement)

Run the code `Main_calculate_displacement.m` to obtain the displacement between the experimental atomic model and the flat atomic model.

### j. Meron and Skyrmion Caculation

Folder: [Calculate_Meron_Skyrmion](./j_Calculate_Meron_Skyrmion)

Run the 6 codes inside the folder `Main_calculate_antimeron_lower.m`, `Main_calculate_antimeron_upper_1.m`, `Main_calculate_antimeron_upper_2.m`, `Main_calculate_meron_lower.m`, `Main_calculate_meron_upper.m`, and `Main_calculate_skyrmion_lower.m` to obatin the merons-like and skyrmion-like structure in the experimental atomic model as well as their corresponding skyrmion number.

### k. Twisted Bilayer Graphene Multislice Simulation with 64 Frozon Phonons

Folder: [Twisted_bilayer_graphene_multislice_simulation_with_64_frozon_phonons](./k_TBG_64_frozen_phonon_simulation)

Run the codes inside the folder `Main_1_AET_RESIRE_Probe.m`, `Main_2_polynomial_tracing_iteration.m`, `Main_3_Match_Atom_With_Position.m`, and `Main_4_Atom_Position_refinement.m` to obatin the traced atom positions for the twisted bilayer graphene multislice simulation with 64 frozon phonons.

### l. Twisted Bilayer Graphene Multislice Simulation with 32 Frozon Phonons

Folder: [Twisted_bilayer_graphene_multislice_simulation_with_32_frozon_phonons](./l_TBG_32_frozen_phonon_simulation)

Run the codes inside the folder `Main_1_AET_RESIRE_Probe.m`, `Main_2_polynomial_tracing_iteration.m`, `Main_3_Match_Atom_With_Position.m`, and `Main_4_Atom_Position_refinement.m` to obatin the traced atom positions for the twisted bilayer graphene multislice simulation with 32 frozon phonons.

### m. Nitrogen Vacancy Center Multislice Simulation

Folder: [Nitrogen_vacancy_center_simulation](./m_Nitrogen_vacancy_center_simulation)

Run the codes inside the folder `Main_1_RESIRE_Multislice_With_Probe.m`, `Main_1_RESIRE_Multislice_Without_Probe.m`, `Main_1_RESIRE_Singleslice_With_Probe.m`, `Main_1_RESIRE_Singleslice_Without_Probe.m`, `Main_2_polynomial_tracing_iteration_Multislice_With_Probe.m`, `Main_2_polynomial_tracing_iteration_Multislice_Without_Probe.m`, `Main_2_polynomial_tracing_iteration_Singleslice_With_Probe.m`, `Main_4_Match_Atom_With_Position_Multislice.m`, `Main_1_RESIRE_Singleslice_With_Probe.m`, `Main_4_Match_Atom_With_Position_Multislice_Without_Probe.m`, `Main_5_Atom_Position_refinement_Multislice_With_Probe.m`, and `Main_5_Atom_Position_refinement_Multislice_Without_Probe.m` to obatin the traced atom positions for the nitrogen vacancy center simulation.

### n. Amorphous Carbon Multislice Simulation

Folder: [Amorphous_carbon_simulation](./n_Amorphous_carbon_simulation)

Run the codes inside the folder `Main_1_AET_RESIRE_Probe.m`,  `Main_2_polynomial_tracing_iteration.m`,  `Main_3_Match_Atom_With_Position.m`, and  `Main_4_Atom_Position_refinement.m` to obatin the traced atom positions for the amorphous carbon multislice simulation.

### o. Resoluction Calculation

Folder: [Calculate_resolution](./o_Calculate_resolution)

Run the codes `Main_calculate_resolution.m` to calculate the resolution of the tomography reconstruction without support constraint.

### p. 2D Tracing Refinement

Folder: [2D_tracing](./p_2D_Tracing)

Run the codes inside the folder `Main_1_AET_RESIRE.m`, `Main_2_2D_Tracing_Lower_Layer.m`, `Main_2_2D_Tracing_Upper_Layer.m`, `Main_3_match_atom_with_projection.m`, and `Main_4_atom_position_refinement.m` to obatin the traced atom positions for the experiment twisted bilayer graphene from 2D Tracing. Run the codes inside the folder `Main_6_match_atom.m`, and `Main_7_calculate_displacement.m` to obtain the corresponding displacement between the experimental atomic model and the flat atomic model.

### p. Odd And Even Dataset Comparison 

Folder: [Experiment_odd_even_dataset](./q_Experiment_odd_even_dataset)

Run the codes inside the Even_dataset and Odd_dataset folders `Main_1_AET_RESIRE_probe.m`, `Main_2_polynomial_tracing_iteration.m`, `Main_3_match_atom_with_projection.m`, `Main_4_atom_position_refinement.m`, `Main_5_match_atom.m`, and `Main_5_Atom_Position_refinement_Multislice_Without_Probe.m` to obatin the traced atom positions for odd and even experiment twisted bilayer graphene dataset and the corresponding displacement between the experimental atomic model and the flat atomic model.



