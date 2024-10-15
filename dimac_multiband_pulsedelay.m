function [pulsedelay,FourierFit1,FourierFit2] = dimac_multiband_pulsedelay(FourierFit1,FourierFit2,R2thresh)
% Function to compare DIMAC pulse waveforms and
% calculate the delay, based on a Fourier series fit to each waveform
% can be child of calc_delay_to_pox.m
%
% Usage:
%       pulsedelay = dimac_pox_pulsedelay(FourierFit1,FourierFit2);
%
%        pulsedelay - output time vector of relative delay for each beat
%      FourierFit1  - object containing the FourierFit output of
%                     dimac_peak_extract.m for the first DIMAC waveform
%      FourierFit2  - object containing the FourierFit output of
%                     dimac_peak_extract.m for the second DIMAC waveform
%               N.B. these two waveforms need to have been acquired
%               simultaneously, so their beat-to-beat timings can be compared
%
%          R2thresh - [optional input - default = 0] R^2 threshold for
%                     processing a pulse, based on fit quality
%
% IDD 29/08/2024
%
% IDD 26/09/2024 - added a check to match number of pulse periods between the two, in cases where they don't match

if nargin < 3
    R2thresh = 0
    % if no threshold specified, does not discard any pulse periods (fit R^2 > 0)
end

%% Calculating the pulse onset for DIMAC (resampled to 1ms time-grid)

% Index each DIMAC timepoint to the scale of a 1ms time-grid:
% (i.e. TR(ms) x interpolated timepoints per DIMAC timepoint)
imgtrigind = 1:FourierFit1.tr*1000:1+(numel(FourierFit1.tc)-1)*FourierFit1.tr*1000;%find(pox_FourierFit.imgtrig(:,1)==1);

FourierFit1.onset_ind = nan(numel(FourierFit1.footind)-1,1);
for beatnum = 1:numel(FourierFit1.footind)-1
    if FourierFit1.R2(beatnum)>R2thresh
        X = fourier_design_matrix(numel(imgtrigind(FourierFit1.footind(beatnum)):imgtrigind(FourierFit1.footind(beatnum+1))),FourierFit1.K,1);
        % Length now based on 1ms sampling period between two DIMAC beats,
        % K is the order of the Fourier coefficients
        beatfit = X*FourierFit1.coeffs(:,beatnum); % Upsampled curve to find the pulse onset
    
        %%%% The following is copied from wrapperPeakFitting.m (fitfoot renamed to onset_ind)
        % calculating fitfoot (the trough preceding the lead edge of the pulse):
        % Window for minimum is the first quarter of the period (between foot indices)
        % Catch for multiple minima with same value in this period is to average
        FourierFit1.onset_ind(beatnum) = imgtrigind(FourierFit1.footind(beatnum)) - 1 + mean(find(beatfit(1:ceil(numel(beatfit)/4))==min(beatfit(1:ceil(numel(beatfit)/4)))));
    end
    
    clear X beatfit
    
end

% Repeat for second DIMAC pulse waveform:

FourierFit2.onset_ind = nan(numel(FourierFit2.footind)-1,1);
for beatnum = 1:numel(FourierFit2.footind)-1
    if FourierFit2.R2(beatnum)>R2thresh
        X = fourier_design_matrix(numel(imgtrigind(FourierFit2.footind(beatnum)):imgtrigind(FourierFit2.footind(beatnum+1))),FourierFit2.K,1);
        % Length now based on 1ms sampling period between two DIMAC beats,
        % K is the order of the Fourier coefficients
        beatfit = X*FourierFit2.coeffs(:,beatnum); % Upsampled curve to find the pulse onset
    
        %%%% The following is copied from wrapperPeakFitting.m (fitfoot renamed to onset_ind)
        % calculating fitfoot (the trough preceding the lead edge of the pulse):
        % Window for minimum is the first quarter of the period (between foot indices)
        % Catch for multiple minima with same value in this period is to average
        FourierFit2.onset_ind(beatnum) = imgtrigind(FourierFit2.footind(beatnum)) - 1 + mean(find(beatfit(1:ceil(numel(beatfit)/4))==min(beatfit(1:ceil(numel(beatfit)/4)))));
    end
    
    clear X beatfit
    
end


if true
    % option to plot onset times on interpolated DIMAC waveforms
    figure
    plot(interp1(imgtrigind,FourierFit1.tc,1:imgtrigind(end)),'k')
    hold on
    plot(interp1(imgtrigind,FourierFit2.tc,1:imgtrigind(end)),'r')
    plot(repmat(FourierFit1.onset_ind',[2 1]),[min(FourierFit1.tc(:))*ones(1,numel(FourierFit1.onset_ind));max(FourierFit1.tc(:))*ones(1,numel(FourierFit1.onset_ind))],'k:')
    plot(repmat(FourierFit2.onset_ind',[2 1]),[min(FourierFit2.tc(:))*ones(1,numel(FourierFit2.onset_ind));max(FourierFit2.tc(:))*ones(1,numel(FourierFit2.onset_ind))],'r:')
    % plot(d_FourierFit.onset_ind(~isnan(d_FourierFit.onset_ind)),pox_FourierFit.imgtrig(d_FourierFit.onset_ind(~isnan(d_FourierFit.onset_ind)),2),'rx')
    % plot(pox_FourierFit.onset_ind(~isnan(pox_FourierFit.onset_ind)),pox_FourierFit.imgtrig(pox_FourierFit.onset_ind(~isnan(pox_FourierFit.onset_ind)),2),'go')
    % plot(pox_FourierFit.footind,pox_FourierFit.imgtrig(pox_FourierFit.footind,2),'b*')
    legend('DIMAC1','DIMAC2','DIMAC1 onset','DIMAC2 onset')
end

%% Calculate the pulse delay from DIMAC1 to DIMAC2:

% IDD 26/09/24 - If different number of points, remove mismatched points by finding the combination of a subset of the larger vector which is closest to the full set of the smaller vector:
if numel(FourierFit1.onset_ind)>numel(FourierFit2.onset_ind)
    crop_ind = nchoosek(1:numel(FourierFit1.onset_ind),numel(FourierFit2.onset_ind));
    FourierFit1.onset_ind = FourierFit1.onset_ind(crop_ind(find(sum((FourierFit1.onset_ind(crop_ind)-repmat(FourierFit2.onset_ind(:)',[size(crop_ind,1) 1])).^2,2,'omitnan')==min(sum((FourierFit1.onset_ind(crop_ind)-repmat(FourierFit2.onset_ind(:)',[size(crop_ind,1) 1])).^2,2,'omitnan')),1),:));
elseif numel(FourierFit1.onset_ind)<numel(FourierFit2.onset_ind)
    crop_ind = nchoosek(1:numel(FourierFit2.onset_ind),numel(FourierFit1.onset_ind));
    FourierFit2.onset_ind = FourierFit2.onset_ind(crop_ind(find(sum((FourierFit2.onset_ind(crop_ind)-repmat(FourierFit1.onset_ind(:)',[size(crop_ind,1) 1])).^2,2,'omitnan')==min(sum((FourierFit2.onset_ind(crop_ind)-repmat(FourierFit1.onset_ind(:)',[size(crop_ind,1) 1])).^2,2,'omitnan')),1),:));
end

pulsedelay = (FourierFit2.onset_ind-FourierFit1.onset_ind)/1000; % Division by 1000 to convert from 1ms indices to delay in seconds

if true
    % option to plot histogram of delay times across beats:
    figure
    histogram(pulsedelay*1000,50)
    xlabel('delay (ms)')
end
