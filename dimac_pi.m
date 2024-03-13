function out = dimac_pi(niiname,out_prefix)
% Usage: out = dimac_pi(niiname,out_prefix);
%
%        niiname  -  full filename of DIMAC file (including extension)
%     out_prefix  -  output prefix name for ROI mask
%
%            out  -  output structure, which will include: DIMAC timeseries (out.ts)
%                                                          Pulsatility Index (out.pi)


addpath(genpath('/home/sapid1/DIMAC_scripts'))

% load in DIMAC dataset
img = load_untouch_nii(niiname);

%% Manual ROI selection:
dimac_process_pipeline(niiname,out_prefix)

mask = load([out_prefix,'.roi.mat']);
% save mask for reference and localisation:
msk1 = img; % Copying header information from DIMAC image
msk1.img = imdilate(mask.roi{1},ones(3,3)); % Replace image data with dilated (4x4) mask
msk1.hdr.dime.dim([1 5]) = [3 1]; % reshape header to 3D
save_untouch_nii(msk1,[out_prefix,'_roi.nii.gz'])



% DIMAC timeseries = mean over 4x4 ROI, centred on manual mask from above (dimac_process_pipeline.m).
ts = mean(reshape(img.img(repmat(imdilate(mask.roi{1},ones(3,3)),[1 1 1 size(img.img,4)])==1),16,size(img.img,4)),1)';
% dilate mask & repmat to 4096 timepoints, so can use to
% index img.img, leaving a 16*4096 vector; reshape
% vector to 16x4096, so can average across the 16
% voxels, giving the mean timecourse

%% Splitting timeseries into heart beat pulse periods:

tr = img.hdr.dime.pixdim(5)
% TR read from the nifti header (units of seconds)

% Peak fitting to split the timeseries into beats:
[tmp.peak,tmp.base,tmp.fit2,tmp.coeffs,tmp.R2,tmp.footind]=dimac_peak_extract(ts,numel(ts),tr);
tmp.tr = tr;
tmp.ts = ts;

%% Calculating pulsatility index (PI) for each beat
% PI = difference between peak and trough of trace, divided by the mean

% Find mean, min and max for each beat (using the raw
% timeseries, rather than the fit, to catch sharp
% peaks)
tmp.beat_mean = zeros(1,numel(tmp.footind)-1);
tmp.beat_max = zeros(1,numel(tmp.footind)-1);
tmp.beat_min = zeros(1,numel(tmp.footind)-1);
for n = 1:numel(tmp.footind)-1
    tmp.beat_mean(n) = mean(ts(tmp.footind(n):tmp.footind(n+1)));
    tmp.beat_max(n) = max(ts(tmp.footind(n):tmp.footind(n+1)));
    tmp.beat_min(n) = min(ts(tmp.footind(n):tmp.footind(n+1)));
end
                    
% Reject beats with poor footind by rejecting beats
% where footind>mean
mean_foot_check = (tmp.fit2(tmp.footind(1:end-1))-tmp.beat_mean')<0;
mean_foot_check(end+1) = (tmp.fit2(tmp.footind(end))-tmp.beat_mean(end))<0; % Checking final footind < mean of preceding beat
mask_mean = mean_foot_check(1:end-1).*mean_foot_check(2:end); % Only include beats where both the start and end footind < mean
                    
% Pulsatility Index
tmp.pi = (tmp.beat_max-tmp.beat_min)./tmp.beat_mean; % Pulsatility Index = range/mean
tmp.pi(mask_mean==0)=nan; % set PI for any rejected beats to NaN
clear mean_foot_check mask_mean

out = tmp;clear tmp


%% Plotting the data, for quality control:

if true
    figure('Position',[100 35 900 485])
    plot(ts)
    hold on
    plot(out.footind,out.fit2(out.footind),'k+')
    plot(out.footind(isnan(out.pi)==0),out.beat_mean(isnan(out.pi)==0),'ms')
    plot(out.footind(isnan(out.pi)==0),out.beat_max(isnan(out.pi)==0),'m^')
    plot(out.footind(isnan(out.pi)==0),out.beat_min(isnan(out.pi)==0),'mv')
    if isnan(mean(out.pi,'omitnan'))==0
        title(['PI (mean/median) = ',num2str(mean(out.pi,'omitnan')),' / ',num2str(median(out.pi,'omitnan'))])
    end
    hold off
end
