function dimac = dimac_pi(dimac)
% IDD 25/09/2024: old version (see 2/2/2024; now dimac_pi_oldROI.m) used old ROI selection method (dimac_process_pipeline)...
%                 ... now updated to accept structure format from dimac_tc.m (13/8/2024) 
%
% Usage: dimac = dimac_pi(dimac);
%
%        dimac    -  structure including dimac timeseries and FourierFit
%        output   -  Pulsatility Index (dimac.pi) added to the structure


addpath(genpath('/home/sapid1/DIMAC_scripts'))

ts = dimac.tc;

%% Calculating pulsatility index (PI) for each beat
% PI = difference between peak and trough of trace, divided by the mean

% Find mean, min and max for each beat (using the raw
% timeseries, rather than the fit, to catch sharp
% peaks)
dimac.beat_mean = zeros(1,numel(dimac.footind)-1);
dimac.beat_max = zeros(1,numel(dimac.footind)-1);
dimac.beat_min = zeros(1,numel(dimac.footind)-1);
for n = 1:numel(dimac.footind)-1
    dimac.beat_mean(n) = mean(ts(dimac.footind(n):dimac.footind(n+1)));
    dimac.beat_max(n) = max(ts(dimac.footind(n):dimac.footind(n+1)));
    dimac.beat_min(n) = min(ts(dimac.footind(n):dimac.footind(n+1)));
end
                    
% Reject beats with poor footind by rejecting beats
% where footind>mean
mean_foot_check = (dimac.fit2(dimac.footind(1:end-1))-dimac.beat_mean')<0;
mean_foot_check(end+1) = (dimac.fit2(dimac.footind(end))-dimac.beat_mean(end))<0; % Checking final footind < mean of preceding beat
mask_mean = mean_foot_check(1:end-1).*mean_foot_check(2:end); % Only include beats where both the start and end footind < mean
                    
% Pulsatility Index
dimac.pi = (dimac.beat_max-dimac.beat_min)./dimac.beat_mean; % Pulsatility Index = range/mean
dimac.pi(mask_mean==0)=nan; % set PI for any rejected beats to NaN
clear mean_foot_check mask_mean



%% Plotting the data, for quality control:

if true
    figure('Position',[100 35 900 485])
    plot(ts)
    hold on
    plot(dimac.footind,dimac.fit2(dimac.footind),'k+')
    plot(dimac.footind(isnan(dimac.pi)==0),dimac.beat_mean(isnan(dimac.pi)==0),'ms')
    plot(dimac.footind(isnan(dimac.pi)==0),dimac.beat_max(isnan(dimac.pi)==0),'m^')
    plot(dimac.footind(isnan(dimac.pi)==0),dimac.beat_min(isnan(dimac.pi)==0),'mv')
    if isnan(mean(dimac.pi,'omitnan'))==0
        title(['PI (mean/median) = ',num2str(mean(dimac.pi,'omitnan')),' / ',num2str(median(dimac.pi,'omitnan'))])
    end
    hold off
end
