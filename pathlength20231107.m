% Notes on calculating the path length between two points on the TOF, based
% on DIMAC ROIs
%
% Test dataset = ISSF10_EFP

% load in data:

addpath(genpath('/home/sapid1/DIMAC_scripts'))

tof = load_untouch_nii('tof_fl3d_tra_8slabs_0p63iso_PAT3.nii.gz');
mca = load_untouch_nii('dpp_Rmca_roi_TOFspace.nii.gz');
car = load_untouch_nii('dpp_Carotid_roi_TOFspace.nii.gz');

% calculate centre of gravity (cog) for each mask:

%carotid:
cog_car(1) = sum(squeeze(sum(sum(car.img,2),3)).*[1:size(car.img,1)]')./sum(squeeze(sum(sum(car.img,2),3)));
cog_car(2) = sum(squeeze(sum(sum(car.img,1),3)).*[1:size(car.img,2)])./sum(squeeze(sum(sum(car.img,1),3)));
cog_car(3) = sum(squeeze(sum(sum(car.img,2),1)).*[1:size(car.img,3)]')./sum(squeeze(sum(sum(car.img,2),1)));
%mca:
cog_mca(1) = sum(squeeze(sum(sum(mca.img,2),3)).*[1:size(mca.img,1)]')./sum(squeeze(sum(sum(mca.img,2),3)));
cog_mca(2) = sum(squeeze(sum(sum(mca.img,1),3)).*[1:size(mca.img,2)])./sum(squeeze(sum(sum(mca.img,1),3)));
cog_mca(3) = sum(squeeze(sum(sum(mca.img,2),1)).*[1:size(mca.img,3)]')./sum(squeeze(sum(sum(mca.img,2),1)));

% Define arteries in the TOF by taking the top 1 percent of voxels as a threshold:
tofmask = tof.img>prctile(tof.img(:),99);

% Set up a starting point at the carotid ROI cog and dilate until it
% overlaps with the artery mask:
se(:,:,1) = [0 0 0;0 1 0;0 0 0]; % se = 6-nearest neighbour structure element 
se(:,:,3) = [0 0 0;0 1 0;0 0 0];
se(:,:,2) = [0 1 0;1 1 1;0 1 0];

startpoint = zeros(size(tofmask));
startpoint(round(cog_car(1)),round(cog_car(2)),round(cog_car(3))) = 1;
while sum(startpoint.*tofmask)==0
    startpoint = imdilate(startpoint,se);
end

% Set up an end point at the MCA cog and dilate until it
% overlaps with the artery mask:
endpoint = zeros(size(tofmask));
endpoint(round(cog_mca(1)),round(cog_mca(2)),round(cog_mca(3))) = 1;
while sum(endpoint.*tofmask)==0
    endpoint = imdilate(endpoint,se);
end

% Crop branches before clustering:
tofmask(:,:,1:round(cog_car(3))-1) = 0; % Crop below the starting point

tofmask(159,170,192) = 0; % Crop Circle of Willis intersection
tofmask(159,169,192) = 0;
tofmask(159,170,189:191) = 0;
tofmask(160,169,189:191) = 0;

tofmask(1:163,:,199:end) = 0; % Crop ACA intersection


% Starting from the carotid (startpoint), dilate amask until it reaches the
% MCA (endpoint). N.B. catch in case it doesn't connect is to break out of
% the loop if the mask size does not increase with an iteration.
amask = startpoint;
s = 0; %Tracking the size of amask on each iteration
while s~=sum(amask(:)) & sum(amask(:).*endpoint(:)) ==0
s = sum(amask(:))
amask = imdilate(amask,se).*tofmask;
end

% Use centre of gravity for each slice to generate a spine through amask:
%%% Did not work, as went to inside of kinks and large error for transverse orientation %%%
% aspine = zeros(size(tof.img));
% for slcnum = round(cog_car(3)):round(cog_mca(3))
%     aspine(round(sum(squeeze(sum(amask(:,:,slcnum),2)).*[1:size(amask,1)]')./sum(squeeze(sum(amask(:,:,slcnum),2)))),round(sum(squeeze(sum(amask(:,:,slcnum),1)).*[1:size(amask,2)])./sum(squeeze(sum(amask(:,:,slcnum),1)))),slcnum) = 1;
% end

% Coordinate matrices, for calculating distances and displacements
x = repmat([1:size(tof.img,1)]',[1 size(tof.img,2) size(tof.img,3)]);
y = repmat([1:size(tof.img,2)],[size(tof.img,1) 1 size(tof.img,3)]);
z = repmat(reshape(1:size(tof.img,3),[1 1 size(tof.img,3)]),[size(tof.img,1) size(tof.img,2) 1]);


stepsize = 2; % distance between search points
initial_vector = [0 0 stepsize]; % Vector to define the next search point
searchdiameter_mm = 5; % 5mm sphere to mask the search for the centre of gravity
searchdiameter_vox = round(mean([searchdiameter_mm./tof.hdr.dime.pixdim(2) searchdiameter_mm./tof.hdr.dime.pixdim(3) searchdiameter_mm./tof.hdr.dime.pixdim(4)]));

% extend amask by the search radius (5mm) past the MCA, so the path can
% extend to the MCA slice position:
endzone = amask;
for count1 = 1:ceil(searchdiameter_vox)
    amask = imdilate(amask,se).*tofmask;
end
endzone = amask-endzone; % Denotes the area where the spine has arrived at the MCA, so the end of the pathlength calculation

aspine = zeros(size(tof.img));
% Initial spine point = centre of gravity at the slice of startpoint:
aspine(round(sum(squeeze(sum(amask(:,:,round(cog_car(3))),2)).*[1:size(amask,1)]')./sum(squeeze(sum(amask(:,:,round(cog_car(3))),2)))),round(sum(squeeze(sum(amask(:,:,round(cog_car(3))),1)).*[1:size(amask,2)])./sum(squeeze(sum(amask(:,:,round(cog_car(3))),1)))),round(cog_car(3))) = 1;

spinevector = [x(aspine==1) y(aspine==1) z(aspine==1)];
next_vector = initial_vector;
count = 1;
while endzone(round(spinevector(count,1)),round(spinevector(count,2)),round(spinevector(count,3)))==0
    
    searchsphere = ((x-(spinevector(count,1)+next_vector(1))).^2+(y-(spinevector(count,2)+next_vector(2))).^2+(z-(spinevector(count,3)+next_vector(3))).^2)<=(searchdiameter_vox/2).^2;
    spinevector(count+1,1) = sum(squeeze(sum(sum(searchsphere.*amask,2),3)).*[1:size(searchsphere.*amask,1)]')./sum(squeeze(sum(sum(searchsphere.*amask,2),3)));
    spinevector(count+1,2) = sum(squeeze(sum(sum(searchsphere.*amask,1),3)).*[1:size(searchsphere.*amask,2)])./sum(squeeze(sum(sum(searchsphere.*amask,1),3)));
    spinevector(count+1,3) = sum(squeeze(sum(sum(searchsphere.*amask,2),1)).*[1:size(searchsphere.*amask,3)]')./sum(squeeze(sum(sum(searchsphere.*amask,2),1)));
    
    next_vector = [spinevector(count+1,1)-spinevector(count,1) spinevector(count+1,2)-spinevector(count,2) spinevector(count+1,3)-spinevector(count,3)];
    next_vector = next_vector.*stepsize./sqrt(sum(next_vector.^2)); % scale next vector so distance = 2 voxels (direction defined by the previous two spinepoints)
    aspine(round(spinevector(count+1,1)),round(spinevector(count+1,2)),round(spinevector(count+1,3))) = 1;
    count = count+1
end

out_pathlength = sum(sqrt(sum(diff(spinevector.*repmat([tof.hdr.dime.pixdim(2) tof.hdr.dime.pixdim(3) tof.hdr.dime.pixdim(4)],size(spinevector,1),1)).^2,2)));
% pathlength taken as the root sum of squares of the difference between
% adjacent points (in mm), summed over all points
