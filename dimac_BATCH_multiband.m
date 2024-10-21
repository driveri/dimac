function [pulse_delay,vessel1,vessel2] = dimac_BATCH_multiband(multiband_nifti_name,tof_nifti_name)
% IDD 10/09/2024: MATLAB function, acts as a batch script (with some BASH) to unpack
%                 a multiband 2-slice nifti DIMAC dataset, define ROIs for two arteries,
%                 extract pulse waveforms for each and calculate the delay between them.
%
%                 N.B.1 this script requires AFNI to be setup, so command-line AFNI commands can be executed
%
% IDD 21/10/2024: Pipeline modified to change how the k-means clustering is performed, so dimac_roi20240809.txt
%                 replaced by dimac_roi20241002.m, which includes a call to dimac_tc.m
%
% Usage:
%
%     [pulse_delay,vessel1,vessel2] = dimac_BATCH_multiband(multiband_nifti_name);
%
%     Inputs:
%         multiband_nifti_name 	- string pointing to the input DIMAC nifti file, including extension
%         tof_nifti_name	- %%%%NOT NEEDED - for future development%%%% string pointing to the TOF nifti file, including extension
%
%     Outputs:
%         pulse_delay	- vector of time delays (in seconds) for each pulse period, between the two arteries
%         vessel1	- [optional] structure including the DIMAC timeseries for the first blood vessel (artery),
%			             the mask used to generate the timeseries and Fourier fitting parameters
%         vessel2	- [optional] equivalent structure for the second blood vessel (artery)

%% Strip file extension:
[filepath1,name1] = fileparts(multiband_nifti_name);
name1 = strtok(name1,'.');

niftiname_no_extension = fullfile(filepath1,name1);
clear filepath1 name1 ext1

%% Split multiband dataset into separate slices
if isfile([multiband_nifti_name,'_slice1.nii.gz'])==0
    eval(['!bash dimac_prepare_multiband.txt ',multiband_nifti_name])
end

%% Runs k-means clustering over pulsatile voxels, then manually select the cluster that corresponds to each artery's pulse waveform and fit the pulse waveform to 5th-order Fourier Series

vessel1 = dimac_roi20241002([niftiname_no_extension,'_slice0.nii.gz'],[niftiname_no_extension,'_slice0_vessel1']);

vessel2 = dimac_roi20241002([niftiname_no_extension,'_slice1.nii.gz'],[niftiname_no_extension,'_slice1_vessel2']);

%% Calculate the delay between the two arteries for each pulse period
[pulse_delay,vessel1,vessel2] = dimac_multiband_pulsedelay(vessel1,vessel2,0.75);
