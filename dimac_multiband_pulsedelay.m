function [pulsedelay,DIMAC1,DIMAC2] = dimac_multiband_pulsedelay(DIMAC1,DIMAC2)
% Function to compare DIMAC pulse waveforms and
% calculate the delay in ms, assuming that onset_ind is in units of the DIMAC TR
%
% Usage:
%       pulsedelay = dimac_pox_pulsedelay(DIMAC1,DIMAC2);
%
%      pulsedelay   - output time vector of relative delay for each beat
%                     (in ms)
%      DIMAC1  - object containing the pulse onset index and TR (in s)
%                     for the first DIMAC waveform
%      DIMAC2  - object containing the pulse onset index and TR (in s)
%                     for the second DIMAC waveform
%               N.B. these two waveforms need to have been acquired
%               simultaneously, so their beat-to-beat timings can be compared
%
%
% IDD 29/08/2024
%
% IDD 26/09/2024 - added a check to match number of pulse periods between the two, in cases where they don't match
% IDD 17/09/2025 - MAJOR CHANGE: adapted script to work for the 2-line fit
%                  onset selection method (Previous version for the Fourier
%                  Series method backed up as
%                  dimac_multiband_pulsedelay_from_FourierSet.m)


if false
    % option to plot onset times on interpolated DIMAC waveforms
    figure
    plot(DIMAC1.tc,'k')
    hold on
    plot(DIMAC2.tc,'r')
    plot(repmat(DIMAC1.onset_ind',[2 1]),[min(DIMAC1.tc(:))*ones(1,numel(DIMAC1.onset_ind));max(DIMAC1.tc(:))*ones(1,numel(DIMAC1.onset_ind))],'k:')
    plot(repmat(DIMAC2.onset_ind',[2 1]),[min(DIMAC2.tc(:))*ones(1,numel(DIMAC2.onset_ind));max(DIMAC2.tc(:))*ones(1,numel(DIMAC2.onset_ind))],'r:')
    legend('DIMAC1','DIMAC2','DIMAC1 onset','DIMAC2 onset')
end

%% Calculate the pulse delay from DIMAC1 to DIMAC2:

% IDD 26/09/24 - If different number of points, remove mismatched points by finding the combination of a subset of the larger vector which is closest to the full set of the smaller vector:
if numel(DIMAC1.onset_ind)>numel(DIMAC2.onset_ind)
    crop_ind = nchoosek(1:numel(DIMAC1.onset_ind),numel(DIMAC2.onset_ind));
    DIMAC1.onset_ind = DIMAC1.onset_ind(crop_ind(find(sum((DIMAC1.onset_ind(crop_ind)-repmat(DIMAC2.onset_ind(:)',[size(crop_ind,1) 1])).^2,2,'omitnan')==min(sum((DIMAC1.onset_ind(crop_ind)-repmat(DIMAC2.onset_ind(:)',[size(crop_ind,1) 1])).^2,2,'omitnan')),1),:));
elseif numel(DIMAC1.onset_ind)<numel(DIMAC2.onset_ind)
    crop_ind = nchoosek(1:numel(DIMAC2.onset_ind),numel(DIMAC1.onset_ind));
    DIMAC2.onset_ind = DIMAC2.onset_ind(crop_ind(find(sum((DIMAC2.onset_ind(crop_ind)-repmat(DIMAC1.onset_ind(:)',[size(crop_ind,1) 1])).^2,2,'omitnan')==min(sum((DIMAC2.onset_ind(crop_ind)-repmat(DIMAC1.onset_ind(:)',[size(crop_ind,1) 1])).^2,2,'omitnan')),1),:));
end

pulsedelay = (DIMAC2.onset_ind-DIMAC1.onset_ind).*DIMAC1.tr*1000; % Conversion from DIMAC TR units to ms delay

if false
    % option to plot histogram of delay times across beats:
    figure
    histogram(pulsedelay,50)
    xlabel('delay (ms)')
end
