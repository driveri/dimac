function [pulsedelay,d_FourierFit,pox_FourierFit] = dimac_pox_pulsedelay(d_FourierFit,pox_FourierFit,R2thresh)
% Function to compare DIMAC and finger pulseox pulse waveforms and
% calculate the delay, based on a Fourier series fit to each waveform
% can be child of calc_delay_to_pox.m
%
% Usage:
%       pulsedelay = dimac_pox_pulsedelay(dimac_FourierFit,pox_FourierFit);
%
%        pulsedelay - output time vector of relative delay for each beat
%      d_FourierFit - object containing the FourierFit output of
%                     dimac_peak_extract.m for the DIMAC waveform
%    pox_FourierFit - object containing the FourierFit output of
%                     dimac_peak_extract.m for the finger pulse waveform
%          R2thresh - [optional input - default = 0] R^2 threshold for
%                     processing a pulse, based on fit quality
%
% IDD 25/01/2024

if nargin < 3
    R2thresh = 0
    % if no threshold specified, does not discard any pulse periods (fit R^2 > 0)
end

%% Calculating the pulse onset for DIMAC (resampled to the pulseox scale)

% Index each DIMAC timepoint to the scale of the pulseox:
imgtrigind = find(pox_FourierFit.imgtrig(:,1)==1);

d_FourierFit.onset_ind = nan(numel(d_FourierFit.footind)-1,1);
for beatnum = 1:numel(d_FourierFit.footind)-1
    if d_FourierFit.R2(beatnum)>R2thresh
        X = fourier_design_matrix(numel(imgtrigind(d_FourierFit.footind(beatnum)):imgtrigind(d_FourierFit.footind(beatnum+1))),d_FourierFit.K,1);
        % Length now based on pulseox sample period between two DIMAC beats, K
        % is the order of the Fourier coefficients
        beatfit = X*d_FourierFit.coeffs(:,beatnum); % Upsampled curve to find the pulse onset
    
        %%%% The following is copied from wrapperPeakFitting.m (fitfoot renamed to onset_ind)
        % calculating fitfoot (the trough preceding the lead edge of the pulse):
        % Window for minimum is the first quarter of the period (between foot indices)
        % Catch for multiple minima with same value in this period is to average
        d_FourierFit.onset_ind(beatnum) = imgtrigind(d_FourierFit.footind(beatnum)) - 1 + mean(find(beatfit(1:ceil(numel(beatfit)/4))==min(beatfit(1:ceil(numel(beatfit)/4)))));
    end
    
    clear X beatfit
    
end


%% Calculating the pulse onset for the pulsox trace (onset_ind)
% If R^2 > R2thresh then calculate onset_ind (index of minimum of Fourier fit)
% Window for minimum is the first quarter of the period (between foot indices)
% Catch for multiple minima with same value in this period is to average

pox_FourierFit.onset_ind = nan(numel(pox_FourierFit.footind)-1,1);
for n = 1:numel(pox_FourierFit.footind)-1
    if pox_FourierFit.R2(n)>R2thresh

        pox_FourierFit.onset_ind(n) = pox_FourierFit.footind(n)-1+mean(find(pox_FourierFit.fit2(pox_FourierFit.footind(n):ceil(pox_FourierFit.footind(n)+(pox_FourierFit.footind(n+1)-pox_FourierFit.footind(n))/4))==min(pox_FourierFit.fit2(pox_FourierFit.footind(n):ceil(pox_FourierFit.footind(n)+(pox_FourierFit.footind(n+1)-pox_FourierFit.footind(n))/4)))));
    end
end

if true
    % option to plot onset times on the raw pulseox trace
    figure
    plot(pox_FourierFit.imgtrig(:,2),'k')
    hold on
    plot(d_FourierFit.onset_ind(~isnan(d_FourierFit.onset_ind)),pox_FourierFit.imgtrig(d_FourierFit.onset_ind(~isnan(d_FourierFit.onset_ind)),2),'rx')
    plot(pox_FourierFit.onset_ind(~isnan(pox_FourierFit.onset_ind)),pox_FourierFit.imgtrig(pox_FourierFit.onset_ind(~isnan(pox_FourierFit.onset_ind)),2),'go')
    plot(pox_FourierFit.footind,pox_FourierFit.imgtrig(pox_FourierFit.footind,2),'b*')
    legend('pulseox','dimac pulse onset','pulseox pulse onset','pulseox foot index')
end

%% Match up dimac and pulseox onset timepoints to calculate the delay:

% debugging: looking for outliers and misalignments of the DIMAC and
% pulseox pulse periods
if false
    nearestpox_diff_poxind = zeros(numel(d_FourierFit.onset_ind),3);
    for count = 1:numel(d_FourierFit.onset_ind)
        nearestpox_diff_poxind(count,1) = find(((pox_FourierFit.onset_ind-d_FourierFit.onset_ind(count)).^2)==min((pox_FourierFit.onset_ind-d_FourierFit.onset_ind(count)).^2));
        nearestpox_diff_poxind(count,2) = d_FourierFit.onset_ind(count)-pox_FourierFit.onset_ind(nearestpox_diff_poxind(count,1));
        nearestpox_diff_poxind(count,3) = pox_FourierFit.onset_ind(nearestpox_diff_poxind(count,1));
    end
    prctile(nearestpox_diff_poxind(:,2),[5 25 50 75 95])
end
%end debugging

% Match DIMAC onset index to nearest pulseox onset index, then calculate
% delay between them (pulsox time - dimac time), so negative delay expected

d_FourierFit.delayfromreference = nan(size(d_FourierFit.onset_ind));
for count = 1:numel(d_FourierFit.onset_ind)
    % find nearest pulse onset point in pulseox trace:
    nearest_pox = find(((pox_FourierFit.onset_ind-d_FourierFit.onset_ind(count)).^2)==min((pox_FourierFit.onset_ind-d_FourierFit.onset_ind(count)).^2));
    
    % calculate time difference from pulseox onset to DIMAC onset (expected
    % to be negative, as intracranial pulse should lead finger pulse)
    tdiff1 = (d_FourierFit.onset_ind(count)-pox_FourierFit.onset_ind(nearest_pox))./pox_FourierFit.fs;
    
    if (nearest_pox+sign(tdiff1) <=numel(pox_FourierFit.onset_ind)) && (nearest_pox+sign(tdiff1)>0)
        % catch to only search of adjacent pulseox within its indices
        if ~isnan(pox_FourierFit.onset_ind(nearest_pox+sign(tdiff1)))
            d_FourierFit.delayfromreference(count) = tdiff1;
            % only record time difference if next closest pulseox (i.e. in
            % opposite direction from the reference point) passes goodness of
            % fit check, so is not NaN.
        end
    end
end

if true
    % option to plot histogram of delay times across beats:
    figure
    histogram(d_FourierFit.delayfromreference*1000,50)
    xlabel('delay (ms)')
end

pulsedelay = d_FourierFit.delayfromreference;
