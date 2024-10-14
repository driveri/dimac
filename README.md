Ian Driver 10/5/2023

A repository of scripts for use in processing DIMAC data (Whittaker 2022, 10.3389/fnins.2021.795749)

Index to follow...

Main Scripts to call:

    dimac_roi20241002.m 	- script to define a DIMAC ROI based on k-means clustering (k=5) of pulsatile voxels
    dimac_tc.m			- GUI for user input into choosing an ROI to generate a DIMAC timeseries
    dimac_pi.m			- script to calculate pulsatililty index, accepting input structure from dimac_tc.m

    readSiemensPhysio.m		- Script to read in Siemens PMU puls data and align with the image timestamps in the DICOM directory
    calc_phys_regressors.m	- Peak fitting script for physiological monitoring traces, written by Kevin Murphy


Dependency scripts (called in the above):

    importSiemens_PULS.m	- Dependency of readSiemensPhysio.m; loads Siemens PMU puls log file

    pulsepowermap.m 		- called in dimac_roi20241002.m; mapping wrapper for pulsepower.mdimac_tc
    pulsepower.m 		- calculates the ratio of power in 40-120 bpm range, scaled by the total power; called by dimac_tc.m and pulsepowermap.m

    overlay_mask.m		- image plotting tool; called in dimac_tc.m GUI
    overlay_jet.m		- image plotting tool; called in dimac_tc.m GUI

    connectedfun2D.m		- 8-nearest neighbour cluster connectivity; called in dimac_tc.m
    dimac_peak_extract.m 	- DIMAC timeseries fitting of Fourier Series (5x sine/cosine pairs)
    fourier_design_matrix.m 	- called by dimac_peak_extract.m and dimac_pox_pulsedelay.m


MRITools scripts for loading and saving nifti (from MRITools 4.0.6):

    load_untouch_nii.m
    load_nii_hdr.m
    load_untouch_nii_hdr.m
    load_nii_ext.m
    load_untouch_nii_img.m
    save_untouch_nii.m
    save_untouch_nii_hdr.m
