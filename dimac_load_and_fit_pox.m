function pox_FourierFit = dimac_load_and_fit_pox(pox,tr)
% Function to load pulseox pulse waveforms.
%
% IDD 28/11/2023
%
% Usage:   pox_FourierFit = dimac_load_and_fit_pox(pox,tr);
%
%      Inputs:   pox             -  string for the pulseox filename (with extension)
%                tr              -  DIMAC TR (in seconds), used to confirm the pulseox sampling rate, based on the scanner trigger spacing
%
%      Output:   pox_FourierFit  -  Structure including pulseox data, scanner trigger positions and peak fitting


%% load pox data and calculate Fourier coefficients:

if isstruct(pox)
    pox_FourierFit = pox;
    
elseif ischar(pox)
    % Load in pulseox trace and scan triggers for carotid scan:
    imgtrig = dlmread(pox);
    
    % find the sample frequency of the phys log, based on the trigger
    % spacing:
    if round(mean(diff(find(imgtrig(:,1))))) == round(tr*1000)
        fs = 1000
    elseif round(2*mean(diff(find(imgtrig(:,1))))) == round(tr*1000)
        fs = 500
    elseif round(2.5*mean(diff(find(imgtrig(:,1))))) == round(tr*1000)
        fs = 400
    else
        error(['Check sample rate and triggers for pox (',num2str(round(mean(diff(find(imgtrig(:,1)))))),').']) 
    end
    
    % Run dimac_peak_extract.m on pulseox trace, then calculate fitfoot:
    [pox_FourierFit.peak,pox_FourierFit.base,pox_FourierFit.fit2,pox_FourierFit.coeffs,pox_FourierFit.R2,pox_FourierFit.footind]=dimac_peak_extract(imgtrig(:,2),numel(imgtrig(:,2)),1/fs);
    % add pulseox trace and scan triggers to structure:
    pox_FourierFit.imgtrig = imgtrig;
    pox_FourierFit.fs = fs; % sample frequency of pulseox trace
    clear fs
    clear imgtrig
    
else
    disp('WARNING: No Pulseox data loaded, check input')
end

pox_FourierFit.K = 5; % Currently hard-coded into dimac_peak_extract.m; K is the # of Fourier pairs fit to each beat.

