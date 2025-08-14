% function to align DIMAC onset indices to pulseox timescales and compare
% DIMAC and pulseox onset times
%
% IDD 11/08/2025
%
% USAGE:    onset_time_delay = calc_onset_delay_dimac_pox(onset_ind_dimac,onset_ind_pox,imgtrig,pox_fs,plotopt);
%
%           onset_time_delay    -   output vector of time delays (in ms) from DIMAC to pulsox pulse onsets (one entry per beat)
%
%           onset_ind_dimac     -   dimac pulse onset index (output of dimac_2linefit_onset.m)
%           onset_ind_pox       -   pulseox pulse onset index (output of dimac_2linefit_onset.m)
%           imgtrig             -   column vector marking where image triggers occurred in respect to the pulseox time indices (1 for trigger, 0 elsewhere)
%           pox_fs              -   pulseox sampling frequency (Hz), used to scale the delay to ms (assumes 1000 Hz sampling if omitted)
%           plotopt             -   logical option to plot a histogram of resulting delay time across beats (only plots if TRUE)

function onset_time_delay = calc_onset_delay_dimac_pox(onset_ind_dimac,onset_ind_pox,imgtrig,pox_fs,plotopt)

trigind = find(imgtrig(:,1));

%% use the scanner trigger index to align the DIMAC pulse onset timings to the pulseox timescale
dimac_oi_in_pox = trigind(floor(onset_ind_dimac))+mod(onset_ind_dimac,1)*diff(trigind(floor(onset_ind_dimac):floor(onset_ind_dimac)+1));


%% For each DIMAC pulse onset: find the closest pulseox pulse onset and calculate the delay from DIMAC to pulsox
onset_time_delay = nan(size(dimac_oi_in_pox));

for n = 1:numel(dimac_oi_in_pox)
    onset_time_delay(n) = onset_ind_pox(abs(onset_ind_pox-dimac_oi_in_pox(n))==min(abs(onset_ind_pox-dimac_oi_in_pox(n))))-dimac_oi_in_pox(n);
end

%% Convert delay from pulseox timepoints to time in ms:
if nargin > 3
    onset_time_delay = onset_time_delay*1000/pox_fs;
end

%% Plot histogram of delay if plotopt TRUE

if nargin > 4
    if plotopt
        figure;histogram(onset_time_delay)
    end
end