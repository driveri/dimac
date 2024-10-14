function out = readSiemensPhysio(PULSfname,DICOMfolder,out_tsv_name)
% Code to read readSiemensPhysio and generate volume triggers, for use with
% peak detection and rebinning of the data based on the cardiac cycle
% N.B. pulse trace interpolated up to 2.5ms temporal resolution, to match
% the Siemens PMU tic resolution, which the scanner triggers are locked to.
%
% IDD 07/10/2024: added an option to output to tsv (based on optional 3rd
%                 input)
%
% Usage:    out = readSiemensPhysio(PULSfname,DICOMfolder,out_tsv_name)
%
% Inputs:       PULSfname      - input Siemens PULS logfile (string)
%               DICOMfolder    - (string) location of the DICOM files, needed
%                                to generate trigger markers for each timepoint
%               out_tsv_name   - (optional; string; no extension) output
%                                filename, this will be two columns, triggers 
%                                and pulse trace, at 2.5 ms resolution
%
% Output:               out    - 3x columns:
%                                       scanner tics (2.5 ms units after midnight)
%                                       scanner triggers, marking when each measurement acquired
%                                       pulse trace

[ACQ_TIME_TICS,CHANNEL,VALUE,SIGNAL] = importSiemens_PULS(PULSfname);
% [VOLUME,SLICE,ACQ_START_TICS,ACQ_FINISH_TICS,ECHO] = importSiemens_Info('Physio_20230113_153045_6d8c998b-f869-4b40-9a06-b7fb6657a3bf_Info.log');

% % Define a volume trigger based on the start time of each TR
% % (ACQ_START_TICS). N.B. where this lies in between two PULS datapoints,
% % this defaults to the later point (PULS sampled at 2 Siemens ticks, so some points lie at the tick in-between)
% vol_trigger = zeros(size(VALUE));
% for count = 1:numel(ACQ_START_TICS)
% vol_trigger(max(find(abs(ACQ_TIME_TICS-ACQ_START_TICS(count))==min(abs(ACQ_TIME_TICS-ACQ_START_TICS(count))))))=1;
% end

% Define a volume trigger based on the DICOM header field 0008,0032 ID Acquisition Time:
dcmsearch = fullfile(DICOMfolder,'*.dcm');
imasearch = fullfile(DICOMfolder,'*.IMA');
dcm_list = [dir(dcmsearch);dir(imasearch)];

acq_time = zeros(numel(dcm_list),1); % Initialise loop vector

for count = 1:numel(dcm_list)
    DICOMfname = dcm_list(count).name;
    dcm_hdr = dicominfo(fullfile(DICOMfolder,DICOMfname));
    acq_time(count) = str2double(dcm_hdr.AcquisitionTime);
end

% Convert from HHMMSS format to PMUtics (units of 2.5ms after midnight):
acq_time_hours = floor(acq_time/10000);
acq_time_minutes = floor((acq_time-floor(acq_time/10000)*10000)/100);
acq_time_seconds = acq_time-floor(acq_time/100)*100;
acq_time_tics = (acq_time_hours*60^2 + acq_time_minutes*60 + acq_time_seconds)*1000/2.5;
% acq_time_tics is equivalent to ACQ_START_TICS from ...INFO.log

% Map volume triggers (acq_time_tics) and pulse waveform (VALUE) onto
% scanner PMU tic scale:

out = zeros(numel(ACQ_TIME_TICS(1):ACQ_TIME_TICS(end)),3);
% 3 columns: [PMUtics,volume_trigger,puls]

out(:,1) = ACQ_TIME_TICS(1):ACQ_TIME_TICS(end);
% PMUtics (N.B. log file records pulse at 2-tic resolution (likely
% interpolated from 20 ms resolution)

out(ismember(out(:,1),round(acq_time_tics)),2)=1;
% Volume trigger = 1; 0 elsewhere

out(:,3) = interp1(ACQ_TIME_TICS,VALUE,out(:,1),'linear');
% Linear interpolation of cardiac pulse waveform (from 2 tics, 5 ms to 1
% tic, 2.5 ms)

% % Option to resample pulse waveform to the TR timescale:
% % N.B. where ACQ_START_TICS lies between two PULS timepoints(ACQ_TIME_TICS)
% % an interpolated signal is calculated by taking the mean of the two
% % adjacent datapoints.
% puls_sync = nan(size(ACQ_START_TICS));
% for count = 1:numel(ACQ_START_TICS)
% puls_sync(count) = mean(VALUE(find(abs(ACQ_TIME_TICS-ACQ_START_TICS(count))==min(abs(ACQ_TIME_TICS-ACQ_START_TICS(count))))));
% end

if nargin > 2
    dlmwrite([out_tsv_name,'.tsv'],out(:,2:3),'\t')
end