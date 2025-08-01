function [onset_ind,fitparams]=dimac_2linefit_onset(x,tr,plotopt)
% Function to take a DIMAC timeseries and calculate the pulse onset time
% based on a 2-line fit of the pulse onset
% 03/06/2025 Ian Driver (IDD)
%
% Usage:
%       [onset_ind,fitparams]=dimac_2linefit_onset(x,tr,plotopt)
%
%       x           - input timeseries
%       tr          - sampling time of the timeseries (in SECONDS)
%       plotopt     - optional input to plot the onset times
%
%       onset_ind   - output vector giving the onset time in ms for each
%                     pulse, based on the first datapoint occuring at 1 ms
%       fitparams   - output structure, capturing the 2-line fit parameters
%                     and R^2 informing the onset time calculation
%                   fitparams.boundind    - start and end boundaries for the
%                                           fitting for each pulse period
%                   fitparams.fit.R2      - R^2 of 2-line fit for each
%                                           vertex point, where max(R^2)
%                                           defines the onset point
%                   fitparams.fit.linfitB - intercept and gradient of
%                                           linear fit up to the vertex
%                                           point
%
% change log:
% DD/MM/YYYY IDD    -   Changed ...

if nargin < 3
    plotopt = false; % default not to make a plot, if option not specified
end

x = double(x); % Catch for load_untouch_nii input, which can be single. This would cause sgolayfilt to error
n = numel(x);

ss = @(x) sum(x.^2); 

%% Peak detection - taken from dimac_peak_extract.m

D = n.*tr;
k = floor(D/3); % sinusoid order for a 3 s cutoff time
X=fourier_design_matrix(n,k,0);
clear k

B=regress(x,X);
x_lp=X*B; % 3s low-pass filter

x_hp=sgolayfilt(x-x_lp,5,21); % Savitzy-Golay filter used to smooth the data before peak detection

% search for peaks bounded by zero crossings after applying a mean+0.75std threshold
mx = -Inf;
startpoint = 1;
endpoint = n;
tempcard = x_hp(startpoint:endpoint);
thresh = mean(tempcard,'omitnan')+0.75*std(tempcard,'omitnan');
%thresh = median(tempcard);
tempcard = ((tempcard-thresh)>0).*(tempcard-thresh);
maxind = [];
k = 1;
for i = 1:length(tempcard)
    if tempcard(i) > mx
        mx = tempcard(i);
        ind = i + startpoint - 1;
    end
    if tempcard(i) == 0 && mx ~= 0;
        maxind(k,1) = ind;
        k = k+1;
        mx = -Inf;
    end
end

% any peaks within 0.5s of a preceding peak are assumed to be local maxima and discarded
Fs=1/tr;
rmThresh=ceil(0.5*Fs);
ii=1;
while (ii < size(maxind,1))
    inddiff=maxind(ii+1,:)-maxind(ii,:);
    while (inddiff < rmThresh)
        maxind(ii+1,:)=[];
        if (ii ~= size(maxind,1))
            inddiff=maxind(ii+1,:)-maxind(ii,:);
        else
            inddiff=Inf;
        end
    end
    ii=ii+1;
end

clear startpoint endpoint tempcard k mx inddiff i ii ind mx X

%% Interpolation up to 1ms:

if round(tr*1e6)-1e3*round(tr*1e3)~=0 % Check that tr is in ms... N.B. rem and mod don't work in this case due to floating point errors
    error('TR is not an integer number of milliseconds, so current interpolation code will fail')
end

t_orig = [1:round(tr/1e-3):round(n*tr/1e-3)]';
t_1ms = [1:round(n*tr/1e-3)]';
x_1ms = interp1(t_orig,x,t_1ms); % linear interpolation
% x_1ms = interpft(x,n*tr/1e-3); % Fourier interpolation

onset_ind = nan(numel(maxind)-1,1);
fitparams.boundind = nan(numel(maxind)-1,2);

%% Loop through pulse periods
for beatnum = 2:numel(maxind)
    
    %% Track back from pulse peak to point of maximum gradient
    % peak position defines the end of the search space
    search_mg_end = find(t_1ms==t_orig(maxind(beatnum)));
    % start of search space is the last quarter of the pulse period
    search_mg_start = find(t_1ms==t_orig(maxind(beatnum)))-ceil((find(t_1ms==t_orig(maxind(beatnum)))-find(t_1ms==t_orig(maxind(beatnum-1))))/4);
    % within this search space, find point of maximum gradient (uses median
    % for case where multiple points have the same max gradient)
    % N.B. as max gradient point is diff between consecutive points, the
    % later is chosen as the reference point maxgradind.
    maxgradind = search_mg_start+floor(median(find(diff(x_1ms(search_mg_start:search_mg_end))==max(diff(x_1ms(search_mg_start:search_mg_end))))));
    clear search_mg_end search_mg_start
    
    %% Start of linear fitting defined by point where curve has similar amplitude to the point of max gradient
    % search space is the first 3/4 of the data
    search_eq_end = find(t_1ms==t_orig(maxind(beatnum)))-ceil((find(t_1ms==t_orig(maxind(beatnum)))-find(t_1ms==t_orig(maxind(beatnum-1))))/4);
    search_eq_start = find(t_1ms==t_orig(maxind(beatnum-1)));
    % find the point with most similar amplitude to that of the point of
    % maximum gradient
    startind = search_eq_start-1+floor(median(find(abs(x_1ms(search_eq_start:search_eq_end)-x_1ms(maxgradind))==min(abs(x_1ms(search_eq_start:search_eq_end)-x_1ms(maxgradind))))));
    clear search_eq_start search_eq_end
    
    fitparams.boundind(beatnum-1,:) = [startind,maxgradind];
    fitparams.fit(beatnum-1).R2 = nan(numel(startind+1:maxgradind-1),1);
    fitparams.fit(beatnum-1).linfitB = nan(numel(startind+1:maxgradind-1),2);
    count1 = 1;
    %% Step through indices, as the boundary (vertex point) between two lines: Optimise boundary as maximum R^2 fit to the data
    for ind1 = startind+1:maxgradind-1
        % linear fit between starting point and boundary
        X = [ones(numel(startind:ind1),1) [startind:ind1]'];
        linfit1 = regress(x_1ms(startind:ind1),X);
        model1 = nan(numel(startind:maxgradind),1);
        model1(1:numel(startind:ind1)) = linfit1(1)+[startind:ind1]*linfit1(2);
        % boundary to end(max gradient point) is a line connecting 
        % the two points
        model1(numel(startind:ind1):end) = linspace(model1(numel(startind:ind1)),x_1ms(maxgradind),numel(ind1:maxgradind));
        tss=ss(x_1ms(startind:maxgradind));
        rss=ss(x_1ms(startind:maxgradind)-model1);
        fitparams.fit(beatnum-1).R2(count1)=1-(rss/tss);% R^2 for the 2-line fit
        
        fitparams.fit(beatnum-1).linfitB(count1,:) = linfit1;
        count1 = count1+1;
        clear X linfit1 model1 tss rss
    end
    onset_ind(beatnum-1) = find(fitparams.fit(beatnum-1).R2==max(fitparams.fit(beatnum-1).R2))+startind;
    clear count1 maxgradind startind
end

if plotopt
    figure('Position',[100 100 800 400])
    plot(t_1ms,x_1ms)
    hold on
    plot(t_orig,x,'k.')
    plot(t_1ms(onset_ind),x_1ms(onset_ind),'rx')
    hold off
    box off
    xlabel('time (ms)')
    ylabel('signal')
end