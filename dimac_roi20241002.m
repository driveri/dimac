function vessel_dimac = dimac_roi20241002(dimac_nifti_name,omask_prefix)
% IDD 02/10/2024: MATLAB function to use k-means clustering (masked by
%                 pulse power ratio) to select an ROI for a vessel
%                 DIMAC pulse waveform
%                 - written to replace dimac_roi20240809.txt, which
%                 restricted the k-means clustering to an area surrounding
%                 the vessel of interest, defined by the TOF
%
% Usage:
%
%     vessel_dimac = dimac_roi20241002(dimac_nifti_name,omask_prefix);
%
%     Inputs:
%         dimac_nifti_name 	- string pointing to the input DIMAC nifti file, including extension
%         omask_prefix	    - [OPTIONAL] string defining the output mask filename
%                             (_roi.nii.gz will be appended to this)
%
%     Output:
%         vessel_dimac      - structure including the DIMAC timeseries for the chosen blood vessel (artery),
%                             the mask used to generate the timeseries and Fourier fitting parameters


%% Strip file extension:
[filepath1,name1] = fileparts(dimac_nifti_name);
name1 = strtok(name1,'.');

niftiname_no_extension = fullfile(filepath1,name1);
clear filepath1 name1 ext1

%% Calculating pulse power ratio (ppr) and restricted k-means clustering:
if isfile([niftiname_no_extension,'_pulsepowerratio.nii.gz'])==0
    pulsepowermap(dimac_nifti_name)
end
if isfile([niftiname_no_extension,'_pulsepowerratio_thresh.nii.gz'])==0
    ppr_threshold = 10; % i.e. 10x as much power in the 40-120bpm frequency range as is in the whole frequency spectrum
    eval(['!3dcalc -a ',niftiname_no_extension,'_pulsepowerratio.nii.gz -expr "step(a-',num2str(ppr_threshold),')" -prefix ',niftiname_no_extension,'_pulsepowerratio_thresh.nii.gz'])
else
    disp(['WARNING: ',niftiname_no_extension,'_pulsepowerratio_thresh.nii.gz already exists, so using original version.'])
end
if isfile([niftiname_no_extension,'_pulsepowerratio_thresh_3dkmeans_k5.nii.gz'])==0
    eval(['!3dkmeans -f ',dimac_nifti_name,' -mask ',niftiname_no_extension,'_pulsepowerratio_thresh.nii.gz -k 5 -prefix ',niftiname_no_extension,'_pulsepowerratio_thresh_3dkmeans_k5.nii.gz'])
else
    disp(['WARNING: ',niftiname_no_extension,'_pulsepowerratio_thresh_3dkmeans_k5.nii.gz already exists, so using original version.'])
end

%% Manual choice of k-means clusters to use and selecting the correct vessel:
if nargin ==1
    % mask nifti not written in this case
    vessel_dimac = dimac_tc(dimac_nifti_name,[niftiname_no_extension,'_pulsepowerratio_thresh_3dkmeans_k5.nii.gz']);
elseif nargin > 1
    % mask nifti written, as filename specified
    vessel_dimac = dimac_tc(dimac_nifti_name,[niftiname_no_extension,'_pulsepowerratio_thresh_3dkmeans_k5.nii.gz'],omask_prefix);
end