% Script to load in ADI data (LabChart), find scanner triggers from a specific scan
% and return the pulse trace data for that scan
% Ian Driver (IDD) 14/04/2023
%%%%
% USAGE: readADI_forDIMACpuls(fname,trig_colnum,puls_colnum,fs)
%
%    fname 		- input filename
%    trig_colnum 	- column number for scanner triggers
%    puls_colnum 	- column number for pulse trace
%    fs 		- sampling frequency (Hz)
%%%%
%
% Trigger formatting and segmentation from Joe Whittaker's script metric_physioproc_spike.m
%
% IDD 21/7/23: changed script to a function to accept a filename (fname), column numbers of scanner triggers (trig_colnum) and pulse (puls_colnum), and sampling rate (fs=500Hz default)
% IDD 1/7/25:  updated padding if logic (line 57-65), so if the file starts < 20 s before the first trigger, the padding only goes back to the first timepoint (avoiding indexing errors)
function readADI_forDIMACpuls(fname,trig_colnum,puls_colnum,fs)

% fname = 'DIMACtest20230331.txt';

% trig_colnum = 4;
% puls_colnum = 3;
if nargin < 4
    fs=1000; % Sampling rate
end

tPad=fs*20;

sdata = dlmread(fname,'\t',9,0); % Load in Physiology logfile data, ignoring 1st 6 rows (header information)

% format trigger column to convert boxcar data (trigger on toggle) to impulse:
trig=[0; abs(diff(double(sdata(:,trig_colnum)>4)))]; % convert to impulse

% pulse waveform:
puls = sdata(:,puls_colnum);

%% Extract scan blocks based on trigger timings:
% Copied straight from Joe's script:
iTrig=find(trig==1);
blockStarts=[iTrig(1); iTrig(([0; diff(iTrig)] > 1*fs))]; % Assumes any gaps > 1 s are breaks between scans
blockEnds=[iTrig(find([0; diff(iTrig)] > 1*fs)-1); iTrig(end)];  % Assumes any gaps > 1 s are breaks between scans
fprintf(1,'\nBlock lengths are...');
for i=1:size(blockStarts,1)
    blockLengths(i,:)=sum(trig(blockStarts(i):blockEnds(i)));
    fprintf(1,'\n\t...%d',blockLengths(i,:));
end
fprintf(1,'\n')

%% Save pulse and trigger data for selected scan:

% trigger block 6 = oblique scan 8
% blockInd = 6; scanname = 'oblique8';
blockInd = input('Enter block number of interest: '); scanname = input('Enter scan name (prefix of output): ','s');

scanStart=blockStarts(blockInd);
scanEnd=blockEnds(blockInd);

vol=zeros(size(trig,1),1);
vol(scanStart:scanEnd,:)=trig(scanStart:scanEnd,:);


% Padding tPad (default 20s) before and after the scan (or to the end of
% the log file if < 20s after the last trigger; and/or from the start of
% the log file if < 20s before the first trigger)
if ((scanEnd+tPad) > size(vol,1)) && ((scanStart-tPad) < 1)
    idxs=1:size(vol,1);
elseif ((scanEnd+tPad) > size(vol,1))
    idxs=(scanStart-tPad):size(vol,1);
elseif (scanStart-tPad) < 1
    idxs=1:(scanEnd+tPad);
else
    idxs=(scanStart-tPad):(scanEnd+tPad);
end
nrows=length(idxs);

physData = zeros(nrows,2); % 2 columns: [vol,puls]

physData(:,1)=vol(idxs,:);
physData(:,2)=puls(idxs,:);

fname_out = [scanname,'_physio.tsv'];
dlmwrite(fname_out,physData,'\t')
