function [peak,base,fit2,coeffs,R2,footind]=dimac_peak_extract(x,n,tr)
% Original function by Joseph Whittaker
% Adapted by Ian Driver (IDD)
%
% Usage:
%       [peak,base,fit2,coeffs,R2,footind]=dimac_peak_extract(x,n,tr);
%
%       x  - input timeseries for peak fitting
%       n  - number of datapoints to be considered (= length of x)
%       tr - sampling time of the timeseries (in SECONDS)
%
% change log:
% 16/08/2023 IDD    -   Changed the high-pass filter cutoff time to
%                       be 3s, depending on input TR and number of
%                       datapoints, giving a duration of timeseries (D=tr*n)
%                       in seconds. Therefore k = D/3 (k is the number of
%                       sine/cosine pairs used in the Fourier Series
%                       filter, with the highest frequency of sinusoid
%                       being k/D Hz)


K=5;

x = double(x); % Catch for load_untouch_nii input, which can be single. This would cause sgolayfilt to error

ss = @(x) sum(x.^2); 

% IDD 16/08/2013: Changed k from 3 to be set by a cut-off frequency of 3s
%X=fourier_design_matrix(n,3,0);
D = n.*tr;
k = floor(D/3); % sinusoid order for a 3 s cutoff time
X=fourier_design_matrix(n,k,0);
clear k
% end of edit IDD 16/08/2023

B=regress(x,X);
x_lp=X*B;

x_hp=sgolayfilt(x-x_lp,5,21);

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

fit=zeros(n,1);
tmptrace=x_hp(1:maxind(1)-1,:);
np=length(tmptrace);
k=round(np/8);
tmpx=fourier_design_matrix(np,k,0);
tmpb=regress(tmptrace,tmpx);
fit(1:maxind(1)-1,:)=tmpx*tmpb;
for i=1:length(maxind)-1
    tmptrace=x_hp(maxind(i):maxind(i+1)-1,:);
    np=length(tmptrace);
    tmpx=fourier_design_matrix(np,3,0);
    tmpb=regress(tmptrace,tmpx);
    fit(maxind(i):maxind(i+1)-1,:)=tmpx*tmpb;
end
tmptrace=x_hp(maxind(end):end,:);
np=length(tmptrace);
k=round(np/8);
tmpx=fourier_design_matrix(np,k,0);
tmpb=regress(tmptrace,tmpx);
fit(maxind(end):end,:)=tmpx*tmpb;
fit=fit+x_lp;

count=1;
if (maxind(1) > 20)
    tmptrace=fit(1:maxind(1)-1,:);
    [~,minIdx]=findpeaks(-1.*tmptrace);
    if length(minIdx)>1 % IDD 17/08/2023: Catch loop added for when no initial minimum detected (i.e. only a slope precedes the first peak)
        footind(count)=minIdx(end);
        count=count+1;
    end
end
for ii=1:size(maxind)-1
    tmptrace=fit(maxind(ii):maxind(ii+1)-1,:);
    [~,minIdx]=findpeaks(-1.*tmptrace);
    if ((length(tmptrace)-minIdx(end)) < 8)
        minIdx(end)=[];
    else
        footind(count)=maxind(ii)+minIdx(end);
        count=count+1;
    end
end


time=linspace(0,n*tr-tr,n);
count=1;

fit2=zeros(n,1);
base=zeros(n,1);
tmptrace=x(1:footind(1),:);
tmpfit=fit(1:footind(1),:);
np=length(tmptrace);
trend=linspace(tmpfit(1),tmpfit(end),np)';
k=round(np/8);
tmpx=fourier_design_matrix(np,k,1);
tmpb=regress((tmptrace-trend),tmpx);
fit2(1:footind(1),:)=tmpx*tmpb;
base(1:footind(1),:)=trend;
[peakVal,peakIdx]=max(tmpx*tmpb);
time1(count,:)=time(peakIdx);
peak1(count,:)=peakVal;
count=count+1;
for i=1:length(footind)-1
    tmptrace=x(footind(i):footind(i+1),:);
    tmpfit=fit(footind(i):footind(i+1),:);
    np=length(tmptrace);
    trend=linspace(tmpfit(1),tmpfit(end),np)';
    tmpx=fourier_design_matrix(np,K,1);
    tmpb=regress((tmptrace-trend),tmpx);
    fit2(footind(i):footind(i+1),:)=tmpx*tmpb;
    base(footind(i):footind(i+1),:)=trend;
    [peakVal,peakIdx]=max(tmpx*tmpb);
    time1(count,:)=time(footind(i)+peakIdx);
    peak1(count,:)=peakVal;
    count=count+1;
    coeffs(:,i)=tmpb(:);
    tss=ss(tmptrace-trend);
    rss=ss((tmptrace-trend)-(tmpx*tmpb));
    R2(:,i)=1-(rss/tss);
end
tmptrace=x(footind(end):end,:);
tmpfit=fit(footind(end):end,:);
np=length(tmptrace);
trend=linspace(tmpfit(1),tmpfit(end),np)';
k=round(np/8);
tmpx=fourier_design_matrix(np,k,1);
tmpb=regress((tmptrace-trend),tmpx);
fit2(footind(end):end,:)=tmpx*tmpb;
base(footind(end):end,:)=trend;
[peakVal,peakIdx]=max(tmpx*tmpb);
time1(count,:)=time(footind(end)+peakIdx);
peak1(count,:)=peakVal;

fit2=fit2+base;

base(1:footind(1),:)=fit(footind(1));
base(footind(end):end,:)=trend(1);
peak1(1)=peak1(2);
peak1(end)=peak1(length(peak1)-1);

peak1=[peak1(1); peak1; peak1(end)];
if time1(1)==0
    time1 = time1(2:end); % Catch to avoid two overlapping timepoints at 0 (which will break the interp call)
    peak1 = peak1(2:end);
end
time1=[0; time1; time(end)];

peak=interp1(time1,peak1,time)';

peak=peak+base;

