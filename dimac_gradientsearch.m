function out = dimac_gradientsearch(tcs,tr)
% Take an array where time is in the first dimension and multiple voxels
% are in higher dimensions

tcs = reshape(tcs,size(tcs,1),[]); % Reshape into a t x #voxels array
tcs = double(tcs); % Catch for load_untouch_nii input, which can be single. This would cause sgolayfilt to error

%% Peak detection - taken from dimac_peak_extract.m and applied to the mean timeseries across voxels:

x = mean(tcs,2);

n = numel(x);

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
    if tempcard(i) == 0 && mx ~= 0
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

clear startpoint endpoint tempcard k mx inddiff i ii ind mx X x

%% Minimum detection - minimum point in the last quarter of each intrabeat period

minind = nan(numel(maxind)-1,size(tcs,2));
gradfit = nan(numel(maxind)-1,size(tcs,2),2);
gradallfit = nan(numel(maxind)-1,size(tcs,2),2);
gradratio = nan(numel(maxind)-1,size(tcs,2));
for voxloop = 1:size(tcs,2)
    x = tcs(:,voxloop);
    for beatnum = 2:numel(maxind)
        search_mg_end = maxind(beatnum);
        search_mg_start = maxind(beatnum)-ceil((maxind(beatnum)-maxind(beatnum-1))/4);
        minind(beatnum-1,voxloop) = search_mg_start+find(x(search_mg_start:search_mg_end)==min(x(search_mg_start:search_mg_end)), 1, 'last' )-1; % Where two points with the same minimum value, choose the later point
        clear search_mg_end search_mg_start
        
        
        %% Calculate gradient from halfway between beats to the minimum point by linear regression
        
        search_grad_start = maxind(beatnum)-ceil((maxind(beatnum)-maxind(beatnum-1))/2);
        search_grad_end = minind(beatnum-1,voxloop);
        
        X1 = [ones(numel(search_grad_start:search_grad_end),1) [search_grad_start:search_grad_end]'];
        X2 = [ones(numel(maxind(beatnum-1):search_grad_end),1) [maxind(beatnum-1):search_grad_end]'];
        
        
        linfit1 = regress(x(search_grad_start:search_grad_end),X1);
        linfit2 = regress(x(maxind(beatnum-1):search_grad_end),X2);
        
        gradfit(beatnum-1,voxloop,:) = linfit1; % Storing the linear fit gradient
        gradallfit(beatnum-1,voxloop,:) = linfit2;
        gradratio(beatnum-1,voxloop) = linfit1(2)./linfit2(2); % Comparing the local gradient to the total gradient, including the preceding peak
        clear search_grad_start search_grad_end X1 X2 linfit1 linfit2
    end
    
    
    
    
end

%% Select voxels where the gradient is negative and the ratio between the two fit gradients is between 0.4 and 2 for over half the beats:

out.mask = sum((gradratio>0.4).*(gradratio<2).*(gradfit(:,:,2)<0),1)'>(numel(maxind)-1)/2;
out.maxind = maxind;
out.minind = minind;
out.gradfit = gradfit;
out.gradallfit = gradallfit;
out.gradratio = gradratio;