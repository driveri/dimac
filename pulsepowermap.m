function pulsepowermap(niftiname,tr)

% Wrapper script for pulsepower.m; loads in the input image and loops
% through voxels, calculating the pulse power ratio for each voxel, then
% saves the output map as a nifti (with _pulsepowerratio appended to the
% input filename)
%
% IDD 04/08/2023
%
% Usage: pulsepowermap(niftiname,tr)
%
%   niftiname   - filename(or complete path) of input DIMAC data (nifti)
%          tr   - optional input for TR (in seconds). If not included, it
%                 is taken from the nifti header.

[niftiname,ext1] = strtok(niftiname,'.'); % Strips the file extension, if included
                               % TO DO:
                               %    - Use fileparts to split the path, in
                               %      case '.' is included in the path structure
                               %    - Use if exist(fname,'file') == 2
                               %      to check that the nifti file exists
                               %      to check for various extensions
                               %      in the case where no extension
                               %      provided in the input
                               
% addpath(genpath('/home/sapid1/DIMAC_scripts'))
nii = load_untouch_nii([niftiname,ext1]);

% If tr optional input not provided, then search in the header for it:
if nargin < 2
    tr = nii.hdr.dime.pixdim(5);
end
disp(['TR = ',num2str(tr),' s, (',num2str(tr*1000),' ms)'])

% Looping over voxels to calculate pulse power ratio:
plspwrratio = zeros(size(nii.img,1),size(nii.img,2),size(nii.img,3));
for x = 1:size(nii.img,1)
    for y = 1:size(nii.img,2)
        for z = 1:size(nii.img,3)
            plspwrratio(x,y,z) = pulsepower(squeeze(double(nii.img(x,y,z,:))),tr);
        end
    end
end

% Preparing nifti structure for saving the output:
nii.img = plspwrratio;
nii.hdr.dime.dim([1 5]) = [3 1];
nii.hdr.dime.datatype=64; % Change the datatype, so not saving output as integers
nii.hdr.dime.bitpix=64;
save_untouch_nii(nii,[niftiname,'_pulsepowerratio',ext1])
