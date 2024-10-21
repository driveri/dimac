Ian Driver 10/5/2023

A repository of scripts for use in processing DIMAC data (Whittaker 2022, 10.3389/fnins.2021.795749)

Index to follow...

Main Scripts to call:

    dimac_roi20241002.m 	- script to define a DIMAC ROI based on k-means clustering (k=5) of pulsatile voxels
    dimac_tc.m			- GUI for user input into choosing an ROI to generate a DIMAC timeseries
    dimac_pi.m			- script to calculate pulsatililty index, accepting input structure from dimac_tc.m

    dimac_prepare_multiband.txt - BASH script to split a 2-slice multiband dataset into individual slices. Use this before running dimac_roi20241002.m on multiband data.
    dimac_BATCH_multiband.m	- script to be used for batch processing - combines the 2-slice multiband DIMAC pipeline to calculate the delay from an inferior slice vessel to a superior slice vessel

    dimac_load_and_fit_pox.m	- loads pulseox log file in format for further processing
    dimac_pox_pulsedelay.m	- calculates the delay between a DIMAC waveform and an associate pulseox (finger pulse) waveform
    dimac_multiband_pulsedelay.m- where 2-slice multiband data available, this skips the pulseox reference and calculates the delay directly between arteries in the 2 slices

    readADI_forDIMACpuls.m 	- Script to read in pulseox data from the CUBRIC ADI physiologicical monitoring setup
    readSiemensPhysio.m		- Script to read in Siemens PMU puls data and align with the image timestamps in the DICOM directory
    calc_phys_regressors.m	- Peak fitting script for physiological monitoring traces, written by Kevin Murphy

    pathlength20231107.m	- Script to mask the TOF image and track the path between two input masks. N.B. this needs to be modified for each TOF, to cut branches


Dependency scripts (called in the above):

    importSiemens_PULS.m	- Dependency of readSiemensPhysio.m; loads Siemens PMU puls log file

    pulsepowermap.m 		- called in dimac_roi20241002.m; mapping wrapper for pulsepower.mdimac_tc
    pulsepower.m 		- calculates the ratio of power in 40-120 bpm range, scaled by the total power; called by dimac_tc.m and pulsepowermap.m

    overlay_mask.m		- image plotting tool; called in dimac_tc.m GUI
    overlay_jet.m		- image plotting tool; called in dimac_tc.m GUI

    connectedfun2D.m		- 8-nearest neighbour cluster connectivity; called in dimac_tc.m
    dimac_peak_extract.m 	- DIMAC timeseries fitting of Fourier Series (5x sine/cosine pairs); called by dimac_tc.m and dimac_load_and_fit_pox.m
    fourier_design_matrix.m 	- called by dimac_peak_extract.m and dimac_pox_pulsedelay.m


MRITools scripts for loading and saving nifti (from MRITools 4.0.6):

    load_untouch_nii.m
    load_nii_hdr.m
    load_untouch_nii_hdr.m
    load_nii_ext.m
    load_untouch_nii_img.m
    save_untouch_nii.m
    save_untouch_nii_hdr.m
