function calc_phys_regressors(filename,prefix,doRetro,doHR,doRVT,doCO2,doO2,retro_ref_file)

%
% Usage:
%    calc_phys_regressors(filename,prefix,doRetro,doHR,doRVT,doCO2,doO2,retro_ref_file)
%				
% 			filename: Name of physiology input file
% 			prefix:   Prefix of the output files
% 			doRetro:  Set to 1 to output RETROICOR regressors 
% 			doHR:     Set to 1 to output heart rate regressors
% 			doRVT:    Set to 1 to output RVT regressors
% 			doCO2:    Set to 1 to output CO2 regressors
% 			doO2:     Set to 1 to output O2 regressors
%		Optional:
% 			retro_ref_file: Name of reference NIFTI file for retroicor output. If this option
% 			                is specified the retroicor regressors will be saved to NIFTI files 
% 			                instead of text files. These outputs can then be used as voxel-
% 			                dependent regressors in an analysis. The purpose of this reference
% 			                file is to provide voxel size, resolution and orientation information.
%
%		 This script will read in physiology files and output regressors to be used in further
%		 analyses. Physiological traces should be stored in columns as text. The program expects
%		 volume and excite triggers to be included also (0 = no trigger, 1 = trigger). Many control 
%		 constants are set within the script. These will be printed each time the program is run. 
%		 Check them carefully and change them if required!!!
%
%		 Written by Kevin Murphy
%		 Last updated 2nd Aug 2013
%		 Email questions to murphyk2@cardiff.ac.uk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Control Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 400; 		% Sampling rate of input physiology file in Hz
ex_col = 2;		% Column in physiology file containing excite triggers (required)
vol_col = 2;	% Column in physiology file containing volume triggers (required)
pox_col = 3;	% Column in physiology file containing pulse oximeter trace (required for doRetro and doHR)
resp_col = 6;	% Column in physiology file containing respiration belt trace (required for doRetro and doRVT)
CO2_col = 9;	% Column in physiology file containing end-tidal CO2 trace (required for doCO2 and doO2)
O2_col = 8;		% Column in physiology file containing end-tidal O2 trace (required for doO2)


respinvert = 1; 			% Set to 1 if the respiration trace is stored inverted 
num_card_terms = 2;		%	Number of RETROICOR cardiac terms to output
num_resp_terms = 2;		%	Number of RETROICOR respiration terms to output
num_int_terms = 1;		%	Number of RETROICOR interaction terms to output (Both num_card_terms and num_resp_terms must be > 0 for this to do anything)
slice_dir = 1; 				%	Indicates slice direction of fMRI prescription: 1 for BottomUp, 2 for TopDown
slice_ord = 1; 				%	Indicates slice order of fMRI prescription: 1 for Interleaved, 2 for sequential
slice_int_first = 1; 	%	The first slice of an interleaved acquistion: either 1 or 2
output_tr = 3;				% TR for the output RETROICOR NIFTI files

CRFconv = 1;	 % Set to 1 to convolve the HR regressor with the CRF and output to separate file
RRFconv = 1;	 % Set to 1 to convolve the RVT regressor with the RRF and output to separate file
HRFconv = 1;	 % Set to 1 to convolve the CO2 and O2 regressors with a standard HRF and output to separate file
output_hr = 1; % Set to 1 to save the full resolution regressors (i.e. at the sampling frequency) along with the convolved
               % regressors and the temporal derivative if requested and the volume triggers; this option is useful for creating different lags

HR_td = 1;	 % Set to 1 to also output the temporal derivative of the HR regressor
RVT_td = 1;	 % Set to 1 to also output the temporal derivative of the RVT regressor
CO2_td = 0;	 % Set to 1 to also output the temporal derivative of the CO2 regressor
O2_td = 1;	 % Set to 1 to also output the temporal derivative of the O2 regressor

demean = 1;    % Set to 1 to demean the HR, RVT, CO2 and O2 regressors and their temporal derivatives before writing them out

ws = 30; % Window size for HR peak detection algorithm (only change this if you know what you are doing!)

peak_check_card = 1; % Set to 1 to turn on plots to check cardiac peak detection accuracy (RECOMMENDED)
peak_check_resp = 1; % Set to 1 to turn on plots to check respiratory peak detection accuracy (RECOMMENDED)
peak_check_CO2 = 1; % Set to 1 to turn on plots to check CO2 peak detection accuracy (RECOMMENDED)
peak_check_O2 = 1; % Set to 1 to turn on plots to check O2 peak detection accuracy (RECOMMENDED)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Control Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Do not change anything below this line %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if (nargin ~= 7 && nargin ~= 8)
   help calc_phys_regressors
   return
end

dos = doRetro + doHR + doRVT + doCO2 + doO2;
if (dos == 0); fprintf(1,'\n\nNo output regressors requested! Exiting.\n\n'); return; end

fprintf(1,'\n\nRunning the calc_phys_regressor script with ');
if (doRetro); fprintf(1,'RETROICOR'); end
if (doRetro && (doHR+doRVT+doCO2>0)); fprintf(1,', '); end
if (doHR); fprintf(1,'HR'); end
if (doHR && (doRVT+doCO2>0)); fprintf(1,', '); end
if (doRVT); fprintf(1,'RVT'); end
if (doRVT && (doCO2>0)); fprintf(1,', '); end
if (doCO2); fprintf(1,'CO2'); end
if (dos > 1); fprintf(1,' and '); end
if (doO2); fprintf(1,'O2'); end
fprintf(1,' regressor outputs.\n');
fprintf(1,'\tInput physiology file: %s\n\tOutput prefix: %s\n',filename,prefix);
if (nargin == 7); fprintf(1,'\n'); end
if (nargin == 8); fprintf(1,'\tRETROICOR reference file: %s\n\n',retro_ref_file); end
fprintf(1,'The following control constants are being used:\n ');
fprintf(1,'\tfs = %d\t\t-\tSampling rate of input physiology file in Hz\n',fs);
fprintf(1,'\tex_col = %d\t\t-\tColumn in physiology file containing excite triggers\n',ex_col);
fprintf(1,'\tvol_col = %d\t\t-\tColumn in physiology file containing volume triggers\n',vol_col);
if (doRetro || doHR)
	fprintf(1,'\tpox_col = %d\t\t-\tColumn in physiology file containing pulse oximeter trace\n',pox_col);
end
if (doRetro || doRVT)
	fprintf(1,'\tresp_col = %d\t\t-\tColumn in physiology file containing respiration belt trace\n',resp_col);
end
if (doRetro || doRVT)
	fprintf(1,'\tresp_col = %d\t\t-\tColumn in physiology file containing respiration belt trace\n',resp_col);
end
if (doCO2 || doO2)
	fprintf(1,'\tCO2_col = %d\t\t-\tColumn in physiology file containing end-tidal CO2 trace\n',CO2_col);
end
if (doO2)
	fprintf(1,'\tO2_col = %d\t\t-\tColumn in physiology file containing end-tidal O2 trace\n',O2_col);
end
fprintf(1,'\n');
if (doRetro || doRVT)
	fprintf(1,'\trespinvert = %d\t\t-\tIf set to 1, the respiration trace will be inverted before processing\n',respinvert);	
end
if (doRetro)
	fprintf(1,'\tnum_card_terms = %d\t-\tNumber of RETROICOR cardiac terms to output\n',num_card_terms);	
	fprintf(1,'\tnum_resp_terms = %d\t-\tNumber of RETROICOR respiration terms to output\n',num_resp_terms);	
	fprintf(1,'\tnum_int_terms = %d\t-\tNumber of RETROICOR interaction terms to output\n',num_int_terms);	
	fprintf(1,'\tslice_dir = %d\t\t-\tIndicates slice direction of fMRI prescription: 1 for BottomUp, 2 for TopDown\n',slice_dir);	
	fprintf(1,'\tslice_ord = %d\t\t-\tIndicates slice order of fMRI prescription: 1 for Interleaved, 2 for sequential\n',slice_ord);	
	if (nargin ==8); fprintf(1,'\toutput_tr = %d\t\t-\tTR for the output RETROICOR NIFTI files\n',output_tr);	end
	fprintf(1,'\n');
end
if (doHR)
	fprintf(1,'\tCRFconv = %d\t\t-\tIf set to 1, the HR regressor convolved with the CRF will also be saved\n',CRFconv);
end
if (doRVT)
	fprintf(1,'\tRRFconv = %d\t\t-\tIf set to 1, the RVT regressor convolved with the RRF will also be saved\n',RRFconv);
end
if (doCO2 || doO2)
	fprintf(1,'\tHRFconv = %d\t\t-\tIf set to 1, the ',HRFconv);
	switch (doCO2 + 2*doO2)
		case 1; fprintf(1,'CO2 regressor ');  
		case 2; fprintf(1,'O2 regressor ');
		case 3; fprintf(1,'CO2 and O2 regressors ');
	end
	fprintf(1,'convolved with the RRF will also be saved\n');
end
if (output_hr)
	fprintf(1,'\toutput_hr = %d\t\t-\tHigh resolution regressors will also be saved\n',output_hr);
end
if (doHR)
	fprintf(1,'\n\tHR_td = %d\t\t-\tIf set to 1, the temporal derivative of the HR regressor will be saved\n',HR_td);
end
if (doRVT)
	fprintf(1,'\tRVT_td = %d\t\t-\tIf set to 1, the temporal derivative of the RVT regressor will be saved\n',RVT_td);
end
if (doCO2)
	fprintf(1,'\tCO2_td = %d\t\t-\tIf set to 1, the temporal derivative of the CO2 regressor will be saved\n',CO2_td);
end
if (doO2)
	fprintf(1,'\tO2_td = %d\t\t-\tIf set to 1, the temporal derivative of the O2 regressor will be saved\n',O2_td);
end
if (demean)
	fprintf(1,'\n\tdemean = %d\t\t-\tIf set to 1, the HR, RVT, CO2 and O2 regressors and their temporal derivatives will be demeaned\n',demean);
end
fprintf(1,'\nRemember to check these constants carefully!!!\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load file and extract relevant columns %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (exist(filename,'file')==0)
	fprintf(1,'\nError! Input physiology file %s does not exist. Exiting\n\n',filename);
	return
end

if (nargin == 8 && exist(retro_ref_file,'file')==0)
	fprintf(1,'\nError! RETROICOR reference file %s does not exist. Exiting\n\n',retro_ref_file);
	return
end


fprintf(1,'Loading physiology file\n');
A = load(filename);
ex = A(:,ex_col);
iexs = find(ex==1);
vol = A(:,vol_col);
ivols = find(vol==1);
if (doRetro || doHR)
	pox = A(:,pox_col);
end
if (doRetro || doRVT)
	resp = A(:,resp_col);
end
if (doO2)
	O2 = A(:,O2_col);
end
if (doCO2 || doO2)
	CO2 = A(:,CO2_col);
end
clear A;

num_vols = length(ivols);
num_sl = length(iexs)/num_vols;

if ((round(num_sl) - num_sl) ~= 0)
	fprintf(1,'\nError:\n\tNumber of excites (%d) not evenly divisable by number of slices (%d)\n\n',length(iexs),length(ivols));
	return
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cardiac Peak Detection Algorithm %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ((doRetro && num_card_terms > 0) || doHR)
	fprintf(1,'Detecting cardiac peaks\n');


	%%%% Old way, not working well when fs~=500Hz
	%	f=[0 0.5 0.7 1  2  3  4 5 6 (fs/2)]/(fs/2);
	%	a=[ 0 0 0.7 1 1 1 0.7 0 0 0];
	%	b=firls(500,f,a);
	%	card=filtfilt(b,1,pox);
	
	frame = round(fs/20); if (round(frame./2) == (frame./2)); frame = frame + 1; end;
	card = sgolayfilt(pox,1,frame);


	maxind = [];
	k = 1;

	for j = 1:floor(length(card)/(fs*ws))
		mx = -Inf;
		startpoint = (j-1)*fs*ws+1;
		endpoint = j*fs*ws;
		tempcard = card(startpoint:endpoint);

		thresh = mean(tempcard)+std(tempcard);
		tempcard = ((tempcard-thresh)>0).*(tempcard-thresh);

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
	end
	
		
	if (peak_check_card)
		%%% Check peak detection results
		x = zeros(length(card),1)+min(card);
		x(maxind) = max(card);
		fig = figure;
  	set(fig,'units','normalized','outerposition',[0 0 1 1]);
		set(fig,'Toolbar','figure');		
		plot(card)
		title('CHANGE CARDIAC PEAKS','FontSize',16)
		hold on
		red = plot(x,'r');
		
%		temp_maxenv = [];
%		temp_maxenv(1:maxind(1)) = card(maxind(1));
%		for j = 1:length(maxind)-1;
%			temp_maxenv(maxind(j):maxind(j+1)) = ...
%				card(maxind(j)) + (card(maxind(j+1))-card(maxind(j))) ...
%						* (((maxind(j):maxind(j+1))-maxind(j)) / (maxind(j+1) - maxind(j)));
%		end
%		temp_maxenv(maxind(end):length(card))=card(maxind(end));
%		black = plot(temp_maxenv,'k');

		but1 = uicontrol('Style','togglebutton','String','Remove Point','units','normalized','Position',[0.01 0.55 0.1 0.05]);
		but2 = uicontrol('Style','togglebutton','String','Insert Point','units','normalized','Position',[0.01 0.5 0.1 0.05]);
		but3 = uicontrol('Style','togglebutton','String','DONE!','units','normalized','Position',[0.01 0.45 0.1 0.05]);

		loop = 1;
		while (loop == 1);
			while ((get(but1,'Value') == 0) && (get(but2,'Value') == 0) && (get(but3,'Value') == 0))
				pause(0.01)
			end
			if (get(but1,'Value') == 1)
				h = impoint;
				pos = getPosition(h);
				delete(h);
				remove_ind = knnsearch(maxind,pos(1));
				maxind(remove_ind) = [];
				x = zeros(length(card),1)+min(card);
				x(maxind) = max(card);
				delete(red);
				red = plot(x,'r');
%				delete(black);
%				temp_maxenv = [];
%				temp_maxenv(1:maxind(1)) = card(maxind(1));
%				for j = 1:length(maxind)-1;
%					temp_maxenv(maxind(j):maxind(j+1)) = ...
%						card(maxind(j)) + (card(maxind(j+1))-card(maxind(j))) ...
%								* (((maxind(j):maxind(j+1))-maxind(j)) / (maxind(j+1) - maxind(j)));
%				end
%				temp_maxenv(maxind(end):length(card))=card(maxind(end));
%				black = plot(temp_maxenv,'k');
				delete(but1);
				but1 = uicontrol('Style','togglebutton','String','Remove Point','units','normalized','Position',[0.01 0.55 0.1 0.05]);
			elseif (get(but2,'Value') == 1)
				h = impoint;
				pos = getPosition(h);
				delete(h);
				insert_ind = round(pos(1));
				maxind = sort([maxind; insert_ind]);
				x = zeros(length(card),1)+min(card);
				x(maxind) = max(card);
				delete(red);
				red = plot(x,'r');
%				delete(black);
%				temp_maxenv = [];
%				temp_maxenv(1:maxind(1)) = card(maxind(1));
%				for j = 1:length(maxind)-1;
%					temp_maxenv(maxind(j):maxind(j+1)) = ...
%						card(maxind(j)) + (card(maxind(j+1))-card(maxind(j))) ...
%								* (((maxind(j):maxind(j+1))-maxind(j)) / (maxind(j+1) - maxind(j)));
%				end
%				temp_maxenv(maxind(end):length(card))=card(maxind(end));
%				black = plot(temp_maxenv,'k');
				delete(but2);
				but2 = uicontrol('Style','togglebutton','String','Insert Point','units','normalized','Position',[0.01 0.5 0.1 0.05]);
			elseif (get(but3,'Value') == 1)
				loop = 0;
				delete(but1);
				delete(but2);
				delete(but3);
				title('Cardiac peaks','FontSize',16)
			end
		end

		pause(0.1)

		output_filename = strcat(prefix,'_card_peakcheck.txt');
		fprintf(1,'\tWriting cardiac peaks to output file: %s\n',output_filename);
		fid=fopen(output_filename,'wt');
		volpeaks = zeros(length(card),1);
		volpeaks(maxind) = 1;
			fprintf(fid,'%f\t%d\n',[card volpeaks]');
	end
	
	
	icardpeaks = maxind;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Respiration Peak Detection Algorithm %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ((doRetro && num_resp_terms > 0) || doRVT)
	fprintf(1,'Detecting respiration peaks .');

	if (respinvert)
		resp = -resp;
	end

	%%%% Old way, not working well when fs~=500Hz
	%  f=[0 0.5/(fs/2) 1/(fs/2)  1];
	%  a=[1 0.98 0 0 ];
  %	 %b=firls(500*5,f,a);
	%  b=firls(500*2,f,a);
	%  respfilt=filtfilt(b,1,resp);
	%  fprintf(1,'.');
	%  
	%  f=[0 0.5/(fs/2) 1/(fs/2)  1];
	%  a=[1 0.6 0 0 ];
	%  %b=firls(500*5,f,a);
	%  b=firls(500*2,f,a);
	%  respfilt1=filtfilt(b,1,resp);
	%  fprintf(1,'.');

	frame = round(fs/2); if (round(frame./2) == (frame./2)); frame = frame + 1; end;
	respfilt = sgolayfilt(resp,1,frame);
	fprintf(1,'.');

	frame = round(fs*2); if (round(frame./2) == (frame./2)); frame = frame + 1; end;
	respfilt1 = sgolayfilt(resp,1,frame);
	fprintf(1,'.');


	%%%%%%%%%%%%%%%
	%%%figure
	%%%plot(resp)
	%%%hold on
	%%%plot(respfilt,'r')
	%%%plot(respfilt1,'g')
	%%%%%%%%%%%%%%%


	diffresp = diff(respfilt1);
			%%% This is a fudge to deal with Anja Hayen's 1000Hz data.
			%%% The problem arises because the respiration trace is not
			%%% truly a 1000Hz signal (I think) and therefore, there is 
			%%% lots of high frequency noise in the diffresp signal. This
			%%% just filters the diffresp again. It's possible that this
			%%% should be done on all data, not just Anja's.
			if (fs == 1000)
				diffresp = sgolayfilt(diffresp,1,frame);
			end
	diffdiffresp = diff(diffresp);

	iminandimax = find(abs(diffresp)<prctile(abs(diffresp),20));
	

	%%% Remove blocks of consecutive indices
	z = zeros(length(respfilt),1);
	z(iminandimax) = 1;
	if (z(1) == 1); z(1) = 0; end;
  inds1 = find(diff(z)==1);
	inds2 = find(diff(z)==-1);
	if (inds2(1) < inds1(1)); inds2(1) = []; end
	iminandimax = [];
	for i = 1:min(length(inds1),length(inds2))
		[m ind] = min(diffresp(inds1(i):inds2(i)));
		iminandimax(i,1) = inds1(i) + ind - 1;
	end
	fprintf(1,'.');

	
	%%%% Remove spurious max and min indices
	%%%i=1;
	%%%while i < length(iminandimax)
	%%%	if iminandimax(i+1)-iminandimax(i) < 0.5*fs
	%%%		if abs(diffresp(iminandimax(i+1))) < abs(diffresp(iminandimax(i)))
	%%%			iminandimax(i) = [];
	%%%		else
	%%%			iminandimax(i+1) = [];
	%%%		end
	%%%		i = i-1;
	%%%	end
	%%%	i = i+1;
	%%%end
	%%%
	%%%if iminandimax(end) >= length(vol)-2
	%%%	iminandimax(end) = [];
	%%%end
	%%%fprintf(1,'.');



	% Calculating maxs
	imax = [];
	imax = find(diffdiffresp(iminandimax)<0);
	imax = iminandimax(imax);
	% Correct for spurious max values
	count = 1;
	while count ~= 0
		maxs = respfilt(imax);
		thresh = mean(maxs) - 2.5*std(maxs);
		icorrect = find(maxs < thresh);
		count = length(icorrect);
		if icorrect ~= 0
			imax(icorrect) = [];
		end
	end

	% Adjust maximums to closest respfilt max
	for i = 1:length(imax)
		dir = 0;
		if (respfilt(imax(i)-1) > respfilt(imax(i))); dir = -1; end;
		if (respfilt(imax(i)+1) > respfilt(imax(i))); dir = 1; end;
		
		while dir ~= 0
			if (respfilt(imax(i)+dir) > respfilt(imax(i)))
				imax(i) = imax(i) + dir;
			else
				dir = 0;
			end
		end
	end
	fprintf(1,'.');


	% Calculating mins
	imin = [];
	for i = 1:length(imax)-1
		[m ind] = min(respfilt(imax(i):imax(i+1)));
		imin(i,1) = imax(i) + ind - 1;
	end
	[m ind] = min(respfilt(imax(end):end));
	imin = [imin; imax(end) + ind - 1];

	fprintf(1,'.\n');

	% Remove maxs and mins that are equal
	iremove = [];
	for i = 1:length(imax)
		if imax(i) == imin(i); iremove = [iremove; i]; end;
	end
	imax(iremove) = [];
	imin(iremove) = [];

	% Remove maxs and mins that within 0.2s of each other
	iremove = [];
	for i = 1:length(imax)
		if abs(imax(i) - imin(i)) < 0.2*fs; iremove = [iremove; i]; end;
	end
	imax(iremove) = [];
	imin(iremove) = [];
	iremove = [];
	for i = 2:length(imax)
		if abs(imax(i) - imin(i-1)) < 0.2*fs; iremove = [iremove; i]; end;
	end
	imax(iremove) = [];
	imin(iremove) = [];



	if (peak_check_resp)
		cont = 0;
		while (cont == 0)

		%%% Check peak detection results for maximums
		x = ones(length(respfilt),1)*max(respfilt)*1.01;
		x(imax) = mean(respfilt);
		fig = figure;
  	set(fig,'units','normalized','outerposition',[0 0 1 1]);
		set(fig,'Toolbar','figure');		
		plot(respfilt)
		title('CHANGE RESPIRATION MAXIMUMS','FontSize',16)
		hold on
		red = plot(x,'r');
		
		temp_maxenv = [];
		temp_maxenv(1:imax(1)) = respfilt(imax(1));
		for j = 1:length(imax)-1;
			temp_maxenv(imax(j):imax(j+1)) = ...
				respfilt(imax(j)) + (respfilt(imax(j+1))-respfilt(imax(j))) ...
						* (((imax(j):imax(j+1))-imax(j)) / (imax(j+1) - imax(j)));
		end
		temp_maxenv(imax(end):length(respfilt))=respfilt(imax(end));
		black = plot(temp_maxenv,'k');

		but1 = uicontrol('Style','togglebutton','String','Remove Point','units','normalized','Position',[0.01 0.65 0.1 0.05]);
		but2 = uicontrol('Style','togglebutton','String','Insert Point','units','normalized','Position',[0.01 0.6 0.1 0.05]);
		but3 = uicontrol('Style','togglebutton','String','Insert Coord','units','normalized','Position',[0.01 0.55 0.1 0.05]);
		but4 = uicontrol('Style','togglebutton','String','DONE!','units','normalized','Position',[0.01 0.45 0.1 0.05]);

		loop = 1;
		while (loop == 1);
			while ((get(but1,'Value') == 0) && (get(but2,'Value') == 0) && (get(but3,'Value') == 0) && (get(but4,'Value') == 0))
				pause(0.01)
			end
			if (get(but1,'Value') == 1)
				h = impoint;
				pos = getPosition(h);
				delete(h);
				remove_ind = knnsearch(imax,pos(1));
				imax(remove_ind) = [];
				x = ones(length(respfilt),1)*max(respfilt)*1.01;
				x(imax) = mean(respfilt);
				delete(red);
				red = plot(x,'r');
				delete(black);
				temp_maxenv = [];
				temp_maxenv(1:imax(1)) = respfilt(imax(1));
				for j = 1:length(imax)-1;
					temp_maxenv(imax(j):imax(j+1)) = ...
						respfilt(imax(j)) + (respfilt(imax(j+1))-respfilt(imax(j))) ...
								* (((imax(j):imax(j+1))-imax(j)) / (imax(j+1) - imax(j)));
				end
				temp_maxenv(imax(end):length(respfilt))=respfilt(imax(end));
				black = plot(temp_maxenv,'k');
				delete(but1);
				but1 = uicontrol('Style','togglebutton','String','Remove Point','units','normalized','Position',[0.01 0.65 0.1 0.05]);
			elseif (get(but2,'Value') == 1)
				h = impoint;
				pos = getPosition(h);
				delete(h);
				insert_ind = round(pos(1));
				[temp_max temp_ind] = max(respfilt(insert_ind-100:insert_ind+100));
				insert_ind = insert_ind-100 + temp_ind;
				imax = sort([imax; insert_ind]);
				x = ones(length(respfilt),1)*max(respfilt)*1.01;
				x(imax) = mean(respfilt);
				delete(red);
				red = plot(x,'r');
				delete(black);
				temp_maxenv = [];
				temp_maxenv(1:imax(1)) = respfilt(imax(1));
				for j = 1:length(imax)-1;
					temp_maxenv(imax(j):imax(j+1)) = ...
						respfilt(imax(j)) + (respfilt(imax(j+1))-respfilt(imax(j))) ...
								* (((imax(j):imax(j+1))-imax(j)) / (imax(j+1) - imax(j)));
				end
				temp_maxenv(imax(end):length(respfilt))=respfilt(imax(end));
				black = plot(temp_maxenv,'k');
				delete(but2);
				but2 = uicontrol('Style','togglebutton','String','Insert Point','units','normalized','Position',[0.01 0.6 0.1 0.05]);
			elseif (get(but3,'Value') == 1)
				h = impoint;
				pos = getPosition(h);
				delete(h);
				insert_ind = round(pos(1));
				imax = sort([imax; insert_ind]);
				respfilt(insert_ind) = pos(2);
				x = ones(length(respfilt),1)*max(respfilt)*1.01;
				x(imax) = mean(respfilt);
				delete(red);
				red = plot(x,'r');
				delete(black);
				temp_maxenv = [];
				temp_maxenv(1:imax(1)) = respfilt(imax(1));
				for j = 1:length(imax)-1;
					temp_maxenv(imax(j):imax(j+1)) = ...
						respfilt(imax(j)) + (respfilt(imax(j+1))-respfilt(imax(j))) ...
								* (((imax(j):imax(j+1))-imax(j)) / (imax(j+1) - imax(j)));
				end
				temp_maxenv(imax(end):length(respfilt))=respfilt(imax(end));
				black = plot(temp_maxenv,'k');
				delete(but3);
				but3 = uicontrol('Style','togglebutton','String','Insert Point','units','normalized','Position',[0.01 0.55 0.1 0.05]);
			elseif (get(but4,'Value') == 1)
				loop = 0;
				delete(but1);
				delete(but2);
				delete(but3);
				delete(but4);
			end
		end

		pause(0.1)

		%%% Check peak detection results for minimums
		x = ones(length(respfilt),1)*min(respfilt)*0.99;
		x(imin) = mean(respfilt);
  	set(fig,'units','normalized','outerposition',[0 0 1 1]);
		set(fig,'Toolbar','figure');		
%		plot(respfilt)
		title('CHANGE RESPIRATION MINIMUMS','FontSize',16)
		hold on
		red = plot(x,'g');
		
		temp_maxenv = [];
		temp_maxenv(1:imin(1)) = respfilt(imin(1));
		for j = 1:length(imin)-1;
			temp_maxenv(imin(j):imin(j+1)) = ...
				respfilt(imin(j)) + (respfilt(imin(j+1))-respfilt(imin(j))) ...
						* (((imin(j):imin(j+1))-imin(j)) / (imin(j+1) - imin(j)));
		end
		temp_maxenv(imin(end):length(respfilt))=respfilt(imin(end));
		black = plot(temp_maxenv,'k');

		but1 = uicontrol('Style','togglebutton','String','Remove Point','units','normalized','Position',[0.01 0.65 0.1 0.05]);
		but2 = uicontrol('Style','togglebutton','String','Insert Point','units','normalized','Position',[0.01 0.6 0.1 0.05]);
		but3 = uicontrol('Style','togglebutton','String','Insert Coord','units','normalized','Position',[0.01 0.55 0.1 0.05]);
		but4 = uicontrol('Style','togglebutton','String','DONE!','units','normalized','Position',[0.01 0.45 0.1 0.05]);

		loop = 1;
		while (loop == 1);
			while ((get(but1,'Value') == 0) && (get(but2,'Value') == 0) && (get(but3,'Value') == 0) && (get(but4,'Value') == 0))
				pause(0.01)
			end
			if (get(but1,'Value') == 1)
				h = impoint;
				pos = getPosition(h);
				delete(h);
				remove_ind = knnsearch(imin,pos(1));
				imin(remove_ind) = [];
				x = ones(length(respfilt),1)*min(respfilt)*0.99;
				x(imin) = mean(respfilt);
				delete(red);
				red = plot(x,'g');
				delete(black);
				temp_maxenv = [];
				temp_maxenv(1:imin(1)) = respfilt(imin(1));
				for j = 1:length(imin)-1;
					temp_maxenv(imin(j):imin(j+1)) = ...
						respfilt(imin(j)) + (respfilt(imin(j+1))-respfilt(imin(j))) ...
								* (((imin(j):imin(j+1))-imin(j)) / (imin(j+1) - imin(j)));
				end
				temp_maxenv(imin(end):length(respfilt))=respfilt(imin(end));
				black = plot(temp_maxenv,'k');
				delete(but1);
				but1 = uicontrol('Style','togglebutton','String','Remove Point','units','normalized','Position',[0.01 0.65 0.1 0.05]);
			elseif (get(but2,'Value') == 1)
				h = impoint;
				pos = getPosition(h);
				delete(h);
				insert_ind = round(pos(1));
				[temp_max temp_ind] = max(respfilt(insert_ind-100:insert_ind+100));
				insert_ind = insert_ind-100 + temp_ind;
				imin = sort([imin; insert_ind]);
				x = ones(length(respfilt),1)*min(respfilt)*0.99;
				x(imin) = mean(respfilt);
				delete(red);
				red = plot(x,'g');
				delete(black);
				temp_maxenv = [];
				temp_maxenv(1:imin(1)) = respfilt(imin(1));
				for j = 1:length(imin)-1;
					temp_maxenv(imin(j):imin(j+1)) = ...
						respfilt(imin(j)) + (respfilt(imin(j+1))-respfilt(imin(j))) ...
								* (((imin(j):imin(j+1))-imin(j)) / (imin(j+1) - imin(j)));
				end
				temp_maxenv(imin(end):length(respfilt))=respfilt(imin(end));
				black = plot(temp_maxenv,'k');
				delete(but2);
				but2 = uicontrol('Style','togglebutton','String','Insert Point','units','normalized','Position',[0.01 0.6 0.1 0.05]);
			elseif (get(but3,'Value') == 1)
				h = impoint;
				pos = getPosition(h);
				delete(h);
				insert_ind = round(pos(1));
				imin = sort([imin; insert_ind]);
				respfilt(insert_ind) = pos(2);
				x = ones(length(respfilt),1)*min(respfilt)*0.99;
				x(imin) = mean(respfilt);
				delete(red);
				red = plot(x,'g');
				delete(black);
				temp_maxenv = [];
				temp_maxenv(1:imin(1)) = respfilt(imin(1));
				for j = 1:length(imin)-1;
					temp_maxenv(imin(j):imin(j+1)) = ...
						respfilt(imin(j)) + (respfilt(imin(j+1))-respfilt(imin(j))) ...
								* (((imin(j):imin(j+1))-imin(j)) / (imin(j+1) - imin(j)));
				end
				temp_maxenv(imin(end):length(respfilt))=respfilt(imin(end));
				black = plot(temp_maxenv,'k');
				delete(but3);
				but3 = uicontrol('Style','togglebutton','String','Insert Point','units','normalized','Position',[0.01 0.55 0.1 0.05]);
			elseif (get(but4,'Value') == 1)
				loop = 0;
				delete(but1);
				delete(but2);
				delete(but3);
				delete(but4);
				title('Respiration maximums and minimums','FontSize',16)
			end
		end

		pause(0.1)

		err = 0;
		last = min([length(imax) length(imin)]);
		
		%%% Check that maxs and mins are alternating
		if (imax(1) < imin(1))
			for i = 1:length(imax) - 1
				if (sum((imin > imax(i)).*(imin < imax(i+1))) ~= 1)
					err = 1;
				end
			end
		end
		if (imin(1) < imax(1))
			for i = 1:length(imin) - 1
				if (sum((imax > imin(i)).*(imax < imin(i+1))) ~= 1)
					err = 1;
				end
			end
		end


		if err == 0
			cont = 1;
		else
			text(length(respfilt)*0.2,mean(respfilt),'Error: maxs and mins must alternate!','FontSize',24,'BackgroundColor','w');
			pause(2)
		end

		
		
		end % while (cont == 0)



		output_filename = strcat(prefix,'_resp_peakcheck.txt');
		fprintf(1,'\tWriting respiration maxs and mins to output file: %s\n',output_filename);
		fid=fopen(output_filename,'wt');
		vol1 = zeros(length(respfilt),1);
		vol1(imax) = 1;
		vol2 = zeros(length(respfilt),1);
		vol2(imin) = 1;
			fprintf(fid,'%f\t%d\t%d\n',[respfilt vol1 vol2]');
	end

	
	irespmaxpeaks = imax(1:min([length(imax) length(imin)]));
	irespminpeaks = imin(1:min([length(imax) length(imin)]));
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CO2 Peak Detection Algorithm %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (doCO2 || doO2)
	fprintf(1,'Detecting end-tidal peaks\n');

	%%%% Old way, not working well when fs~=500Hz
	%	f=[0 5/fs 10/fs  1];
	%	a=[1 0.87 0 0 ];
	%	b=firls(100,f,a);
	%	CO2filt=filtfilt(b,1,CO2);
	%	O2filt=filtfilt(b,1,O2);

	frame = round(fs/5); if (round(frame./2) == (frame./2)); frame = frame + 1; end;
	CO2filt = sgolayfilt(CO2,1,frame);
	if (doO2); O2filt = sgolayfilt(O2,1,frame); end

	%%%%%%%%%%%%%%%
	%%%figure
	%%%plot(CO2)
	%%%hold on
	%%%plot(CO2filt,'r')
	%%%figure
	%%%plot(O2)
	%%%hold on
	%%%plot(O2filt,'r')
	%%%%%%%%%%%%%%%





	% Calculate where we are going to assume the subject starts breathing out
	x = ((diff(CO2filt)-max(diff(CO2filt))/4)>0).*(diff(CO2filt)-max(diff(CO2filt))/4);
	indzero = find(x==0);
	indrise = find(x(indzero(1:end-1)+1)>0);
	imaxstart = indzero(indrise)+1;
	
	% Calc position of max between the start of one positive lobe and the next
	imax = [];
	[m ind] = max(CO2filt(1:imaxstart(2)));
	imax(1) = ind;
	for j=2:length(imaxstart)-1;
		[m ind] = max(CO2filt(imaxstart(j):imaxstart(j+1)));
		imax(j) = imaxstart(j) + ind;
	end
	[m ind] = max(CO2filt(imaxstart(end):length(CO2filt)));
	imax(length(imaxstart)) = imaxstart(end) + ind;
	imax = imax'-1;


	% Remove spurious values where the max is the same as the value at imaxstart
	iremove = find(imaxstart(2:end)==imax(1:end-1));
	imax(iremove) = [];

	% Remove short breaths where end-tidal CO2 is not reached with a sliding window technique
	iremove = [];
	if (imax(1) == 0); imax(1) = 1; end
	for i = 6:length(imax)-5;
		m = mean([CO2filt(imax(i-5:i-1));CO2filt(imax(i+1:i+5))]);
		l = min(CO2filt(imax(i-5):imax(i+5)));
		if (CO2filt(imax(i)) < (0.5 * (m + l))); iremove = [iremove;i]; end
	end
	imax(iremove) = [];

% Check for peaks that are within 2s of each other and choose only the highest peak
	iremove = [];
	i = 1;
	while (i < length(imax));
		count = 0;
		j = 1;
		check_inds = [];
		too_long = 0;
		while ((too_long == 0) && (imax(i+j) - imax(i+j-1)) < 2 * fs)
			check_inds = [check_inds;i+j-1];
			j = j+1;
			count = 1;
			if ((i+j) > length(imax)); too_long = 1; end
		end
	
		if (count == 1)
			check_inds = [check_inds;i+j-1];
			[m ind] = max(CO2filt(imax(check_inds)));
			check_inds(ind) = [];
			iremove = check_inds;
			imax(iremove) = [];
			imaxstart(iremove) = [];	
		else
			i = i+1; 
		end
	end

	if (peak_check_CO2 && doCO2)
		%%% Check peak detection results
		x = zeros(length(CO2filt),1);
		x(imax) = max(CO2filt);
		fig = figure;
  	set(fig,'units','normalized','outerposition',[0 0 1 1]);
		set(fig,'Toolbar','figure');		
		plot(CO2filt)
		title('CHANGE CO2 PEAKS','FontSize',16)
		hold on
		red = plot(x,'r');
		
		temp_maxenv = [];
		temp_maxenv(1:imax(1)) = CO2filt(imax(1));
		for j = 1:length(imax)-1;
			temp_maxenv(imax(j):imax(j+1)) = ...
				CO2filt(imax(j)) + (CO2filt(imax(j+1))-CO2filt(imax(j))) ...
						* (((imax(j):imax(j+1))-imax(j)) / (imax(j+1) - imax(j)));
		end
		temp_maxenv(imax(end):length(CO2filt))=CO2filt(imax(end));
		black = plot(temp_maxenv,'k');

		but1 = uicontrol('Style','togglebutton','String','Remove Point','units','normalized','Position',[0.01 0.65 0.1 0.05]);
		but2 = uicontrol('Style','togglebutton','String','Insert Point','units','normalized','Position',[0.01 0.6 0.1 0.05]);
		but3 = uicontrol('Style','togglebutton','String','Insert Coord','units','normalized','Position',[0.01 0.55 0.1 0.05]);
		but4 = uicontrol('Style','togglebutton','String','DONE!','units','normalized','Position',[0.01 0.45 0.1 0.05]);
		
		loop = 1;
		while (loop == 1);
			while ((get(but1,'Value') == 0) && (get(but2,'Value') == 0) && (get(but3,'Value') == 0) && (get(but4,'Value') == 0))
				pause(0.01)
			end
			if (get(but1,'Value') == 1)
				h = impoint;
				pos = getPosition(h);
				delete(h);
				remove_ind = knnsearch(imax,pos(1));
				imax(remove_ind) = [];
				x = zeros(length(CO2filt),1);
				x(imax) = max(CO2filt);
				delete(red);
				red = plot(x,'r');
				delete(black);
				temp_maxenv = [];
				temp_maxenv(1:imax(1)) = CO2filt(imax(1));
				for j = 1:length(imax)-1;
					temp_maxenv(imax(j):imax(j+1)) = ...
						CO2filt(imax(j)) + (CO2filt(imax(j+1))-CO2filt(imax(j))) ...
								* (((imax(j):imax(j+1))-imax(j)) / (imax(j+1) - imax(j)));
				end
				temp_maxenv(imax(end):length(CO2filt))=CO2filt(imax(end));
				black = plot(temp_maxenv,'k');
				delete(but1);
				but1 = uicontrol('Style','togglebutton','String','Remove Point','units','normalized','Position',[0.01 0.65 0.1 0.05]);
			elseif (get(but2,'Value') == 1)
				h = impoint;
				pos = getPosition(h);
				delete(h);
				insert_ind = round(pos(1));
				[temp_max temp_ind] = max(CO2filt(insert_ind-100:insert_ind+100));
				insert_ind = insert_ind-100 + temp_ind;
				imax = sort([imax; insert_ind]);
				x = zeros(length(CO2filt),1);
				x(imax) = max(CO2filt);
				delete(red);
				red = plot(x,'r');
				delete(black);
				temp_maxenv = [];
				temp_maxenv(1:imax(1)) = CO2filt(imax(1));
				for j = 1:length(imax)-1;
					temp_maxenv(imax(j):imax(j+1)) = ...
						CO2filt(imax(j)) + (CO2filt(imax(j+1))-CO2filt(imax(j))) ...
								* (((imax(j):imax(j+1))-imax(j)) / (imax(j+1) - imax(j)));
				end
				temp_maxenv(imax(end):length(CO2filt))=CO2filt(imax(end));
				black = plot(temp_maxenv,'k');
				delete(but2);
				but2 = uicontrol('Style','togglebutton','String','Insert Point','units','normalized','Position',[0.01 0.6 0.1 0.05]);
			elseif (get(but3,'Value') == 1)
				h = impoint;
				pos = getPosition(h);
				delete(h);
				insert_ind = round(pos(1));
				imax = sort([imax; insert_ind]);
				CO2filt(insert_ind) = pos(2);
				x = zeros(length(CO2filt),1);
				x(imax) = max(CO2filt);
				delete(red);
				red = plot(x,'r');
				delete(black);
				temp_maxenv = [];
				temp_maxenv(1:imax(1)) = CO2filt(imax(1));
				for j = 1:length(imax)-1;
					temp_maxenv(imax(j):imax(j+1)) = ...
						CO2filt(imax(j)) + (CO2filt(imax(j+1))-CO2filt(imax(j))) ...
								* (((imax(j):imax(j+1))-imax(j)) / (imax(j+1) - imax(j)));
				end
				temp_maxenv(imax(end):length(CO2filt))=CO2filt(imax(end));
				black = plot(temp_maxenv,'k');
				delete(but3);
				but3 = uicontrol('Style','togglebutton','String','Insert Coord','units','normalized','Position',[0.01 0.55 0.1 0.05]);
			elseif (get(but4,'Value') == 1)
				loop = 0;
				delete(but1);
				delete(but2);
				delete(but3);
				delete(but4);
				title('End-tidal CO2 peaks','FontSize',16)
			end
		end

		pause(0.1)
		
		output_filename = strcat(prefix,'_CO2_peakcheck.txt');
		fprintf(1,'\tWriting CO2 peaks to output file: %s\n',output_filename);
		fid=fopen(output_filename,'wt');
		volpeaks = zeros(length(CO2filt),1);
		volpeaks(imax) = 1;
			fprintf(fid,'%f\t%d\n',[CO2filt volpeaks]');
	end

	

	iCO2peaks = imax;



	if (doO2)
		iO2peaks = iCO2peaks;
	
		%%%%%% Move O2 peaks to closest minimum within 1 second
		for i = 1:length(iO2peaks)-1
			[m,ind] = min(O2filt(iO2peaks(i):iO2peaks(i)+fs));
			iO2peaks(i) = iO2peaks(i) + ind;
		end
		if (iO2peaks(end)+fs < length(O2filt))
			[m,ind] = min(O2filt(iO2peaks(end):iO2peaks(end)+fs));
			iO2peaks(end) = iO2peaks(end) + ind;
		end
	end




	if (peak_check_O2 && doO2)
		%%% Check peak detection results
		x = zeros(length(O2filt),1)+min(O2filt);
		x(iO2peaks) = mean(O2filt);
		fig = figure;
  	set(fig,'units','normalized','outerposition',[0 0 1 1]);
		set(fig,'Toolbar','figure');		
		plot(O2filt)
		title('CHANGE O2 TROUGHS','FontSize',16)
		hold on
		red = plot(x,'r');
		
		temp_maxenv = [];
		temp_maxenv(1:iO2peaks(1)) = O2filt(iO2peaks(1));
		for j = 1:length(iO2peaks)-1;
			temp_maxenv(iO2peaks(j):iO2peaks(j+1)) = ...
				O2filt(iO2peaks(j)) + (O2filt(iO2peaks(j+1))-O2filt(iO2peaks(j))) ...
						* (((iO2peaks(j):iO2peaks(j+1))-iO2peaks(j)) / (iO2peaks(j+1) - iO2peaks(j)));
		end
		temp_maxenv(iO2peaks(end):length(O2filt))=O2filt(iO2peaks(end));
		black = plot(temp_maxenv,'k');

		but1 = uicontrol('Style','togglebutton','String','Remove Point','units','normalized','Position',[0.01 0.65 0.1 0.05]);
		but2 = uicontrol('Style','togglebutton','String','Insert Point','units','normalized','Position',[0.01 0.6 0.1 0.05]);
		but3 = uicontrol('Style','togglebutton','String','Insert Coord','units','normalized','Position',[0.01 0.55 0.1 0.05]);
		but4 = uicontrol('Style','togglebutton','String','DONE!','units','normalized','Position',[0.01 0.45 0.1 0.05]);

		loop = 1;
		while (loop == 1);
			while ((get(but1,'Value') == 0) && (get(but2,'Value') == 0) && (get(but3,'Value') == 0) && (get(but4,'Value') == 0))
				pause(0.01)
			end
			if (get(but1,'Value') == 1)
				h = impoint;
				pos = getPosition(h);
				delete(h);
				remove_ind = knnsearch(iO2peaks,pos(1));
				iO2peaks(remove_ind) = [];
				x = zeros(length(O2filt),1)+min(O2filt);
				x(iO2peaks) = mean(O2filt);
				delete(red);
				red = plot(x,'r');
				delete(black);
				temp_maxenv = [];
				temp_maxenv(1:iO2peaks(1)) = O2filt(iO2peaks(1));
				for j = 1:length(iO2peaks)-1;
					temp_maxenv(iO2peaks(j):iO2peaks(j+1)) = ...
						O2filt(iO2peaks(j)) + (O2filt(iO2peaks(j+1))-O2filt(iO2peaks(j))) ...
								* (((iO2peaks(j):iO2peaks(j+1))-iO2peaks(j)) / (iO2peaks(j+1) - iO2peaks(j)));
				end
				temp_maxenv(iO2peaks(end):length(O2filt))=O2filt(iO2peaks(end));
				black = plot(temp_maxenv,'k');
				delete(but1);
				but1 = uicontrol('Style','togglebutton','String','Remove Point','units','normalized','Position',[0.01 0.65 0.1 0.05]);
			elseif (get(but2,'Value') == 1)
				h = impoint;
				pos = getPosition(h);
				delete(h);
				insert_ind = round(pos(1));
				[temp_max temp_ind] = max(O2filt(insert_ind-100:insert_ind+100));
				insert_ind = insert_ind-100 + temp_ind;
				iO2peaks = sort([iO2peaks; insert_ind]);
				x = zeros(length(O2filt),1)+min(O2filt);
				x(iO2peaks) = mean(O2filt);
				delete(red);
				red = plot(x,'r');
				delete(black);
				temp_maxenv = [];
				temp_maxenv(1:iO2peaks(1)) = O2filt(iO2peaks(1));
				for j = 1:length(iO2peaks)-1;
					temp_maxenv(iO2peaks(j):iO2peaks(j+1)) = ...
						O2filt(iO2peaks(j)) + (O2filt(iO2peaks(j+1))-O2filt(iO2peaks(j))) ...
								* (((iO2peaks(j):iO2peaks(j+1))-iO2peaks(j)) / (iO2peaks(j+1) - iO2peaks(j)));
				end
				temp_maxenv(iO2peaks(end):length(O2filt))=O2filt(iO2peaks(end));
				black = plot(temp_maxenv,'k');
				delete(but2);
				but2 = uicontrol('Style','togglebutton','String','Insert Point','units','normalized','Position',[0.01 0.6 0.1 0.05]);
			elseif (get(but3,'Value') == 1)
				h = impoint;
				pos = getPosition(h);
				delete(h);
				insert_ind = round(pos(1));
				iO2peaks = sort([iO2peaks; insert_ind]);
				O2filt(insert_ind) = pos(2);
				x = zeros(length(O2filt),1)+min(O2filt);
				x(iO2peaks) = mean(O2filt);
				delete(red);
				red = plot(x,'r');
				delete(black);
				temp_maxenv = [];
				temp_maxenv(1:iO2peaks(1)) = O2filt(iO2peaks(1));
				for j = 1:length(iO2peaks)-1;
					temp_maxenv(iO2peaks(j):iO2peaks(j+1)) = ...
						O2filt(iO2peaks(j)) + (O2filt(iO2peaks(j+1))-O2filt(iO2peaks(j))) ...
								* (((iO2peaks(j):iO2peaks(j+1))-iO2peaks(j)) / (iO2peaks(j+1) - iO2peaks(j)));
				end
				temp_maxenv(iO2peaks(end):length(O2filt))=O2filt(iO2peaks(end));
				black = plot(temp_maxenv,'k');
				delete(but3);
				but3 = uicontrol('Style','togglebutton','String','Insert Point','units','normalized','Position',[0.01 0.55 0.1 0.05]);
			elseif (get(but4,'Value') == 1)
				loop = 0;
				delete(but1);
				delete(but2);
				delete(but3);
				delete(but4);
				title('End-tidal O2 troughs','FontSize',16)
			end
		end

    pause(0.1)


		output_filename = strcat(prefix,'_O2_peakcheck.txt');
		fprintf(1,'\tWriting O2 peaks to output file: %s\n',output_filename);
		fid=fopen(output_filename,'wt');
		volpeaks = zeros(length(O2filt),1);
		volpeaks(iO2peaks) = 1;
			fprintf(fid,'%f\t%d\n',[O2filt volpeaks]');
	end

	

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate and write out RETROICOR regressors %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (doRetro)
	fprintf(1,'Calculating RETROICOR regressors\n');

	if (exist('icardpeaks') == 0)
		fprintf(1,'\n\nOutput of cardiac peak detection algorithm does not exist\nExiting!\n\n');
		return
	end
	if (exist('respfilt') == 0)
		fprintf(1,'\n\nOutput of respiration filtering does not exist\nExiting!\n\n');
		return
	end
	if (exist('irespmaxpeaks') == 0 || exist('irespminpeaks') == 0)
		fprintf(1,'\n\nOutput of respiration peak detection algorithm does not exist\nExiting!\n\n');
		return
	end

	% Reorder excite triggers so that they are in order of BottomToTop and Sequential
	reorder_iexs = iexs;
	if (slice_dir == 1 && slice_ord == 1 && slice_int_first == 1)
		x = reshape(reorder_iexs,num_sl,num_vols);
		i1 = (1:ceil(num_sl/2))';
		i2 = (ceil(num_sl/2)+1:num_sl)'; 
		reord = [];
		count = 1;
		for j = 1:length(i2)
			reord(count) = i1(j);
			count = count + 1;
			reord(count) = i2(j);
			count = count + 1;
		end
		if (mod(num_sl,2)==1)
			reord(count) = i1(end);
		end
		y = x(reord',:);
		reorder_iexs = y(:);
	elseif (slice_dir == 1 && slice_ord == 1 && slice_int_first == 2)
		x = reshape(reorder_iexs,num_sl,num_vols);
		i1 = (ceil((num_sl+1)/2):num_sl)';
		i2 = (1:ceil((num_sl+1)/2)-1)';
		reord = [];
		count = 1;
		for j = 1:length(i2)
			reord(count) = i1(j);
			count = count + 1;
			reord(count) = i2(j);
			count = count + 1;
		end
		if (mod(num_sl,2)==1)
			reord(count) = i1(end);
		end
		y = x(reord',:);
		reorder_iexs = y(:);
	elseif (slice_dir == 1 && slice_ord == 2)
		% Do Nothing
	elseif (slice_dir == 2 && slice_ord == 1 && slice_int_first == 1)
		x = reshape(reorder_iexs,num_sl,num_vols);
		i1 = (ceil(num_sl/2):-1:1)';
		i2 = (num_sl:-1:ceil(num_sl/2)+1)'; 
		reord = [];
		count = 1;
		if (mod(num_sl,2)==1)
			for j = 1:length(i2)
				reord(count) = i1(j);
				count = count + 1;
				reord(count) = i2(j);
				count = count + 1;
			end
			reord(count) = i1(end);
		else
			for j = 1:length(i2)
				reord(count) = i2(j);
				count = count + 1;
				reord(count) = i1(j);
				count = count + 1;
			end
		end
		y = x(reord',:);
		reorder_iexs = y(:);
	elseif (slice_dir == 2 && slice_ord == 1 && slice_int_first == 2)
		x = reshape(reorder_iexs,num_sl,num_vols);
		i1 = (num_sl:-1:ceil((num_sl-1)/2)+1)';
		i2 = (ceil((num_sl-1)/2):-1:1)';
		reord = [];
		count = 1;
		if (mod(num_sl,2)==1)
			for j = 1:length(i2)
				reord(count) = i1(j);
				count = count + 1;
				reord(count) = i2(j);
				count = count + 1;
			end
			reord(count) = i1(end);
		else
			for j = 1:length(i2)
				reord(count) = i2(j);
				count = count + 1;
				reord(count) = i1(j);
				count = count + 1;
			end
		end
		y = x(reord',:);
		reorder_iexs = y(:);
	elseif (slice_dir == 2 && slice_ord == 2)
		x = reshape(reorder_iexs,num_sl,num_vols);
		reord = (num_sl:-1:1);
		y = x(reord',:);
		reorder_iexs = y(:);
	end





	% Calculate thetac
	if (num_card_terms > 0)
		thetac = zeros(1,length(card));
		j=1;
		startindex=0;
		endindex=icardpeaks(j);
		for i=1:length(card)
 			if i == endindex
   			if j == length(icardpeaks)
     			startindex=icardpeaks(end);
     			endindex=length(card);
   			else
     			startindex=icardpeaks(j);
     			j=j+1;
     			endindex=icardpeaks(j);
   			end
 			end
			thetac(i) = 2 * pi * (i - startindex) / (endindex - startindex);
		end
		cardphases = thetac(reorder_iexs);
		cardphases = reshape(cardphases,length(iexs)/length(ivols),length(ivols))';
	end



	% Calculate thetar
	if (num_resp_terms > 0)
		thetar = zeros(length(iexs),1);
		resptrace = (respfilt(reorder_iexs) - min(respfilt(reorder_iexs)))/abs(max(respfilt(reorder_iexs) - min(respfilt(reorder_iexs))));
		resphist = hist(resptrace,100);
		respdenom = sum(resphist);
	
		respsign = zeros(length(respfilt),1);
		if (irespmaxpeaks(1) < irespminpeaks(1))
			respsign(1:irespmaxpeaks(1)) = 1;
			for i = 1:length(irespmaxpeaks)-1 
				respsign(irespmaxpeaks(i):irespminpeaks(i)) = -1;
				respsign(irespminpeaks(i):irespmaxpeaks(i+1)) = 1;
			end
			respsign(irespmaxpeaks(end):end) = -1;
		else
			respsign(1:irespminpeaks(1)) = -1;
			for i = 1:length(irespmaxpeaks)-1 
				respsign(irespminpeaks(i):irespmaxpeaks(i)) = 1;
				respsign(irespmaxpeaks(i):irespminpeaks(i+1)) = -1;
			end
			respsign(irespminpeaks(end):end) = 1;
		end
		respsign = respsign(reorder_iexs);
	
		for i = 1:length(thetar)
			thetar(i) = pi * sum(resphist(1:round(resptrace(i)*100))) / respdenom .* respsign(i);
		end
		respphases = thetar;
		respphases = reshape(respphases,length(iexs)/length(ivols),length(ivols))';
	end
	

	% Calculate sine and cosine terms
	if (num_card_terms > 0)
		for t = 1:num_card_terms
			CardSinTerms(:,:,t) = sin(t*cardphases);
			CardCosTerms(:,:,t) = cos(t*cardphases);
		end
	end
	if (num_resp_terms > 0)
		for t = 1:num_resp_terms
			RespSinTerms(:,:,t) = sin(t*respphases);
			RespCosTerms(:,:,t) = cos(t*respphases);
		end
	end
	if (num_int_terms > 0 && num_card_terms > 0 && num_resp_terms > 0)
 		for t = 1:num_int_terms
  		for tt = 1:num_int_terms
    		IntSinAddTerm(:,:,t,tt) = sin(t*cardphases+tt*respphases);
     		IntCosAddTerm(:,:,t,tt) = cos(t*cardphases+tt*respphases);
     		IntSinSubTerm(:,:,t,tt) = sin(t*cardphases-tt*respphases);
     		IntCosSubTerm(:,:,t,tt) = cos(t*cardphases-tt*respphases);
  		end
		end
	end


	% Write out RETROICOR text files
	if (nargin ==7)
		if (num_card_terms > 0)
			for t = 1:num_card_terms
				output_filename = strcat(prefix,'_RET_CardSin',num2str(t),'.txt');
				fprintf(1,'\tWriting RETROICOR regressor to output file: %s\n',output_filename);
				x = CardSinTerms(:,:,t);
				save(output_filename,'x','-ASCII','-tabs')
			
				output_filename = strcat(prefix,'_RET_CardCos',num2str(t),'.txt');
				fprintf(1,'\tWriting RETROICOR regressor to output file: %s\n',output_filename);
				x = CardCosTerms(:,:,t);
				save(output_filename,'x','-ASCII','-tabs')
			end
		end

		if (num_resp_terms > 0)
			for t = 1:num_resp_terms
				output_filename = strcat(prefix,'_RET_RespSin',num2str(t),'.txt');
				fprintf(1,'\tWriting RETROICOR regressor to output file: %s\n',output_filename);
				x = RespSinTerms(:,:,t);
				save(output_filename,'x','-ASCII','-tabs')
			
				output_filename = strcat(prefix,'_RET_RespCos',num2str(t),'.txt');
				fprintf(1,'\tWriting RETROICOR regressor to output file: %s\n',output_filename);
				x = RespCosTerms(:,:,t);
				save(output_filename,'x','-ASCII','-tabs')
			end
		end

		if (num_int_terms > 0 && num_card_terms > 0 && num_resp_terms > 0)
 			for t = 1:num_int_terms
  			for tt = 1:num_int_terms
 					output_filename = strcat(prefix,'_RET_IntSin+',num2str(t),num2str(tt),'.txt');
					fprintf(1,'\tWriting RETROICOR regressor to output file: %s\n',output_filename);
					x = IntSinAddTerm(:,:,t);
					save(output_filename,'x','-ASCII','-tabs')

 					output_filename = strcat(prefix,'_RET_IntCos+',num2str(t),num2str(tt),'.txt');
					fprintf(1,'\tWriting RETROICOR regressor to output file: %s\n',output_filename);
					x = IntCosAddTerm(:,:,t);
					save(output_filename,'x','-ASCII','-tabs')

 					output_filename = strcat(prefix,'_RET_IntSin-',num2str(t),num2str(tt),'.txt');
					fprintf(1,'\tWriting RETROICOR regressor to output file: %s\n',output_filename);
					x = IntSinSubTerm(:,:,t);
					save(output_filename,'x','-ASCII','-tabs')

 					output_filename = strcat(prefix,'_RET_IntCos-',num2str(t),num2str(tt),'.txt');
					fprintf(1,'\tWriting RETROICOR regressor to output file: %s\n',output_filename);
					x = IntCosSubTerm(:,:,t);
					save(output_filename,'x','-ASCII','-tabs')

  			end
			end
		end
	end
	
	
	% Write out RETROICOR NIFTI files
	gzip_ref = 0;
	if (nargin == 8)
		fprintf(1,'Reading RETROICOR reference file %s\n',retro_ref_file);
		if (strfind(retro_ref_file,'.gz') == (length(retro_ref_file) - 2))
			[status,result] = system(['gunzip ' retro_ref_file]);
			if (status~=0) 
				fprintf(1,'\nError running gunzip on %s. Exiting!\n\n',retro_ref_file);
				return
			end
			retro_ref_file = retro_ref_file(1:end-3);
			gzip_ref = 1;
		end
		ref_nii = load_untouch_nii(retro_ref_file);
		if (gzip_ref); system(['gzip ' retro_ref_file]); end

		
		new_nii = ref_nii;
		if (find(new_nii.hdr.dime.dim == num_sl) == 4)
			xdim = new_nii.hdr.dime.dim(2);
			ydim = new_nii.hdr.dime.dim(3);
			slice_dim = 3;
		elseif (find(new_nii.hdr.dime.dim == num_sl) == 3)
			xdim = new_nii.hdr.dime.dim(2);
			ydim = new_nii.hdr.dime.dim(4);
			slice_dim = 2;
		elseif (find(new_nii.hdr.dime.dim == num_sl) == 2)
			xdim = new_nii.hdr.dime.dim(3);
			ydim = new_nii.hdr.dime.dim(4);
			slice_dim = 1;
		else
			fprintf(1,'\nError: number of slices in %s not the same as in %s\n',retro_ref_file,filename);
			return
		end
		
		precision = 'double';
		new_nii.hdr.dime.datatype = 64;
		new_nii.hdr.dime.dim(1) = 4;
		new_nii.hdr.dime.dim(5) = num_vols;
		new_nii.hdr.dime.pixdim(5) = output_tr;

		new_nii.img = zeros(xdim,ydim,num_sl,num_vols,precision);
	

		if (num_card_terms > 0)
			for t = 1:num_card_terms
				
				output_filename = strcat(prefix,'_RET_CardSin',num2str(t),'.nii');
				new_nii.fileprefix = strcat(prefix,'_RET_CardSin',num2str(t));
				fprintf(1,'\tWriting RETROICOR regressor to output file: %s.gz\n',output_filename);
				new_nii.img = zeros(xdim,ydim,num_sl,num_vols,precision);
				x = CardSinTerms(:,:,t);
				for i = 1:xdim
					for j = 1:ydim
						new_nii.img(i,j,:,:) = x';
					end
				end
				if (slice_dim == 1); new_nii.img = permute(new_nii.img,[3 1 2 4]); end;
				if (slice_dim == 2); new_nii.img = permute(new_nii.img,[1 3 2 4]); end;
				save_untouch_nii(new_nii,output_filename);
				system(['gzip ' output_filename]);

				output_filename = strcat(prefix,'_RET_CardCos',num2str(t),'.nii');
				new_nii.fileprefix = strcat(prefix,'_RET_CardCos',num2str(t));
				fprintf(1,'\tWriting RETROICOR regressor to output file: %s.gz\n',output_filename);
				new_nii.img = zeros(xdim,ydim,num_sl,num_vols,precision);
				x = CardCosTerms(:,:,t);
				for i = 1:xdim
					for j = 1:ydim
						new_nii.img(i,j,:,:) = x';
					end
				end
				if (slice_dim == 1); new_nii.img = permute(new_nii.img,[3 1 2 4]); end;
				if (slice_dim == 2); new_nii.img = permute(new_nii.img,[1 3 2 4]); end;
				save_untouch_nii(new_nii,output_filename);
				system(['gzip ' output_filename]);
			
			end
		end
		
		
		if (num_resp_terms > 0)
			for t = 1:num_resp_terms
				
				output_filename = strcat(prefix,'_RET_RespSin',num2str(t),'.nii');
				new_nii.fileprefix = strcat(prefix,'_RET_RespSin',num2str(t));
				fprintf(1,'\tWriting RETROICOR regressor to output file: %s.gz\n',output_filename);
				new_nii.img = zeros(xdim,ydim,num_sl,num_vols,precision);
				x = RespSinTerms(:,:,t);
				for i = 1:xdim
					for j = 1:ydim
						new_nii.img(i,j,:,:) = x';
					end
				end
				if (slice_dim == 1); new_nii.img = permute(new_nii.img,[3 1 2 4]); end;
				if (slice_dim == 2); new_nii.img = permute(new_nii.img,[1 3 2 4]); end;
				save_untouch_nii(new_nii,output_filename);
				system(['gzip ' output_filename]);

				output_filename = strcat(prefix,'_RET_RespCos',num2str(t),'.nii');
				new_nii.fileprefix = strcat(prefix,'_RET_RespCos',num2str(t));
				fprintf(1,'\tWriting RETROICOR regressor to output file: %s.gz\n',output_filename);
				new_nii.img = zeros(xdim,ydim,num_sl,num_vols,precision);
				x = RespCosTerms(:,:,t);
				for i = 1:xdim
					for j = 1:ydim
						new_nii.img(i,j,:,:) = x';
					end
				end
				if (slice_dim == 1); new_nii.img = permute(new_nii.img,[3 1 2 4]); end;
				if (slice_dim == 2); new_nii.img = permute(new_nii.img,[1 3 2 4]); end;
				save_untouch_nii(new_nii,output_filename);
				system(['gzip ' output_filename]);
			
			end
		end

		if (num_int_terms > 0)
 			for t = 1:num_int_terms
  			for tt = 1:num_int_terms
				
					output_filename = strcat(prefix,'_RET_IntSin+',num2str(t),num2str(tt),'.nii');
					new_nii.fileprefix = strcat(prefix,'_RET_IntSin+',num2str(t),num2str(tt));
					fprintf(1,'\tWriting RETROICOR regressor to output file: %s.gz\n',output_filename);
					new_nii.img = zeros(xdim,ydim,num_sl,num_vols,precision);
					x = IntSinAddTerm(:,:,t);
					for i = 1:xdim
						for j = 1:ydim
							new_nii.img(i,j,:,:) = x';
						end
					end
					if (slice_dim == 1); new_nii.img = permute(new_nii.img,[3 1 2 4]); end;
					if (slice_dim == 2); new_nii.img = permute(new_nii.img,[1 3 2 4]); end;
					save_untouch_nii(new_nii,output_filename);
					system(['gzip ' output_filename]);

					output_filename = strcat(prefix,'_RET_IntCos+',num2str(t),num2str(tt),'.nii');
					new_nii.fileprefix = strcat(prefix,'_RET_IntCos+',num2str(t),num2str(tt));
					fprintf(1,'\tWriting RETROICOR regressor to output file: %s.gz\n',output_filename);
					new_nii.img = zeros(xdim,ydim,num_sl,num_vols,precision);
					x = IntCosAddTerm(:,:,t);
					for i = 1:xdim
						for j = 1:ydim
							new_nii.img(i,j,:,:) = x';
						end
					end
					if (slice_dim == 1); new_nii.img = permute(new_nii.img,[3 1 2 4]); end;
					if (slice_dim == 2); new_nii.img = permute(new_nii.img,[1 3 2 4]); end;
					save_untouch_nii(new_nii,output_filename);
					system(['gzip ' output_filename]);

					output_filename = strcat(prefix,'_RET_IntSin-',num2str(t),num2str(tt),'.nii');
					new_nii.fileprefix = strcat(prefix,'_RET_IntSin-',num2str(t),num2str(tt));
					fprintf(1,'\tWriting RETROICOR regressor to output file: %s.gz\n',output_filename);
					new_nii.img = zeros(xdim,ydim,num_sl,num_vols,precision);
					x = IntSinSubTerm(:,:,t);
					for i = 1:xdim
						for j = 1:ydim
							new_nii.img(i,j,:,:) = x';
						end
					end
					save_untouch_nii(new_nii,output_filename);
					system(['gzip ' output_filename]);

					output_filename = strcat(prefix,'_RET_IntCos-',num2str(t),num2str(tt),'.nii');
					new_nii.fileprefix = strcat(prefix,'_RET_IntCos-',num2str(t),num2str(tt));
					fprintf(1,'\tWriting RETROICOR regressor to output file: %s.gz\n',output_filename);
					new_nii.img = zeros(xdim,ydim,num_sl,num_vols,precision);
					x = IntCosSubTerm(:,:,t);
					for i = 1:xdim
						for j = 1:ydim
							new_nii.img(i,j,:,:) = x';
						end
					end
					if (slice_dim == 1); new_nii.img = permute(new_nii.img,[3 1 2 4]); end;
					if (slice_dim == 2); new_nii.img = permute(new_nii.img,[1 3 2 4]); end;
					save_untouch_nii(new_nii,output_filename);
					system(['gzip ' output_filename]);

  			end
			end
		end
		
	end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate and write out HR regressor %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (doHR)
	fprintf(1,'Calculating HR regressor\n');

	if (exist('icardpeaks') == 0)
		fprintf(1,'\n\nOutput of cardiac peak detection algorithm does not exist\nExiting!\n\n');
		return
	end

	cardrate = 60./diff(icardpeaks/fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%% REMOVED IN FAVOUR OF USING peak_check_card
%%%
%%%	% Correct for spurious cardiac rates
%%%	count = 1;
%%%	while count ~=0
%%%		count = 0;
%%%		posthresh = mean(cardrate)+3*std(cardrate);
%%%		negthresh = mean(cardrate)-3*std(cardrate);
%%%		ipos = find(cardrate > posthresh);
%%%		ineg = find(cardrate < negthresh);
%%%		if length(ipos) ~= 0
%%%			if ipos(1) == 1
%%%				cardrate(ipos(1)) = cardrate(ipos(1)+1);
%%%				ipos(1) = [];
%%%			end
%%%			if (ipos(end) == length(cardrate))
%%%				cardrate(ipos(end)) = cardrate(ipos(end)-1);
%%%				ipos(end) = [];
%%%			end
%%%			cardrate(ipos) = (cardrate(ipos-1) + cardrate(ipos+1))/2;
%%%			count = 1;
%%%		end
%%%		if length(ineg) ~= 0
%%%			if (ineg(1) == 1)
%%%				cardrate(ineg(1)) = cardrate(ineg(1)+1);
%%%				ineg(1) = [];
%%%			end
%%%			if (length(ineg) ~= 0 && ineg(end) == length(cardrate))
%%%				cardrate(ineg(end)) = cardrate(ineg(end)-1);
%%%				ineg(end) = [];
%%%			end
%%%			cardrate(ineg) = (cardrate(ineg-1) + cardrate(ineg+1))/2;
%%%			count = 1;
%%%		end
%%%	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	HR = zeros(length(card),1);
	HR(1:icardpeaks(1)) = cardrate(1);
	for j=1:length(icardpeaks)-2
		HR(icardpeaks(j):icardpeaks(j+1)) = cardrate(j) + (cardrate(j+1)-cardrate(j)) ...
						* (((icardpeaks(j):icardpeaks(j+1))-icardpeaks(j)) / (icardpeaks(j+1) - icardpeaks(j)));
	end 
	HR(icardpeaks(end-1):length(card))=cardrate(end);

	
	HRregressor = HR(ivols);
	% Make some attempt at smoothing
	HRregressor(2:end-1) = (0.5*HRregressor(1:end-2) + HRregressor(2:end-1) + 0.5*HRregressor(3:end))/2;

	if (HR_td); HRregressor_td = diff(HR); HRregressor_td = HRregressor_td(ivols); end


	output_filename = strcat(prefix,'_HR.txt');
	fprintf(1,'\tWriting HR regressor to output file: %s\n',output_filename);
	fid=fopen(output_filename,'wt');
	if (HR_td)
		HRregressor_td = diff(HR); HRregressor_td = HRregressor_td(ivols);
		if (demean); HRregressor = HRregressor - mean(HRregressor); HRregressor_td = HRregressor_td - mean(HRregressor_td); end
			fprintf(fid,'%f\t%f\n',[HRregressor HRregressor_td]');
	else
		if (demean); HRregressor = HRregressor - mean(HRregressor); end
			fprintf(fid,'%f\n',HRregressor');	
	end
	fclose(fid);

	if (CRFconv == 1)
		fprintf(1,'\tConvolving HR regressor with CRF\n');
		t = 0:1/fs:25;
		CRF = 0.6 .* (t .^ 2.7) .* (exp(-t./1.6)) - 16 .* (1./sqrt(18*pi)) .* exp(-0.5 .* (1/9) .* (t-12) .* (t-12));
		HRconv = conv([HR(1)*ones(60*fs,1);HR],CRF)/sum(CRF);
		HRconv = HRconv(60*fs+1:60*fs+length(HR));
		HRconvregressor = HRconv(ivols);
		if (HR_td); HRconvregressor_td = diff(HRconv); HRconvregressor_td = HRconvregressor_td(ivols); end

		output_filename = strcat(prefix,'_HR_CRFconv.txt');
		fprintf(1,'\tWriting HR CRFconv regressor to output file: %s\n',output_filename);
		fid=fopen(output_filename,'wt');
		if (HR_td)
			if (demean); HRconvregressor = HRconvregressor - mean(HRconvregressor); HRconvregressor_td = HRconvregressor_td - mean(HRconvregressor_td); end
				fprintf(fid,'%f\t%f\n',[HRconvregressor HRconvregressor_td]');
		else
			if (demean); HRconvregressor = HRconvregressor - mean(HRconvregressor); end
				fprintf(fid,'%f\n',HRconvregressor');
		end
		fclose(fid);
	end

	if (output_hr)
		output_filename = strcat(prefix,'_HR_hr.txt');
		fprintf(1,'\tWriting high resolution HR regressor to output file: %s\n',output_filename);
		fid=fopen(output_filename,'wt');
		if (HR_td)
			if (CRFconv)
				HRtd = diff(HRconv); HRtd(end+1) = HRtd(end);
				if (demean); HR = HR - mean(HR); HRconv = HRconv - mean(HRconv(1:length(HR))); HRtd = HRtd - mean(HRtd(1:length(HR))); end
					fprintf(fid,'%f\t%f\t%f\t%d\n',[HR HRconv HRtd vol]');
			else
				HRtd = diff(HR); HRtd(end+1) = HRtd(end);
				if (demean); HR = HR - mean(HR); HRtd = HRtd - mean(HRtd); end
					fprintf(fid,'%f\t%f\t%d\n',[HR HRtd vol]');
			end
		else
			if (CRFconv)
				if (demean); HR = HR - mean(HR); HRconv = HRconv - mean(HRconv(1:length(HR)));end
						fprintf(fid,'%f\t%f\t%d\n',[HR HRconv vol]');
			else
				if (demean); HR = HR - mean(HR); end
					fprintf(fid,'%f\t%d\n',[HR vol]');
			end
		end
		fclose(fid);
	end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate and write out RVT regressor %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (doRVT)
	fprintf(1,'Calculating RVT regressor\n');

	if (exist('irespmaxpeaks') == 0 || exist('irespminpeaks') == 0)
		fprintf(1,'\n\nOutput of respiration peak detection algorithm does not exist\nExiting!\n\n');
		return
	end

	maxtrace = zeros(length(respfilt),1);
	mintrace = zeros(length(respfilt),1);
	maxtrace(1:irespmaxpeaks(1)) = respfilt(irespmaxpeaks(1));
	mintrace(1:irespminpeaks(1)) = respfilt(irespminpeaks(1));
	maxtrace(irespmaxpeaks(end):end) = respfilt(irespmaxpeaks(end));
	mintrace(irespminpeaks(end):end) = respfilt(irespminpeaks(end));

	for j = 1:length(irespmaxpeaks)-1
		maxtrace(irespmaxpeaks(j):irespmaxpeaks(j+1)) = respfilt(irespmaxpeaks(j)) ...
					+ (respfilt(irespmaxpeaks(j+1))-respfilt(irespmaxpeaks(j))) ...
								* (((irespmaxpeaks(j):irespmaxpeaks(j+1))-irespmaxpeaks(j)) / (irespmaxpeaks(j+1) - irespmaxpeaks(j)));
	end

	for j = 1:length(irespminpeaks)-1
		mintrace(irespminpeaks(j):irespminpeaks(j+1)) = respfilt(irespminpeaks(j)) ...
					+ (respfilt(irespminpeaks(j+1))-respfilt(irespminpeaks(j))) ...
								* (((irespminpeaks(j):irespminpeaks(j+1))-irespminpeaks(j)) / (irespminpeaks(j+1) - irespminpeaks(j)));
	end
	
	periodtrace = zeros(length(respfilt),1);
	for j = 1:length(irespmaxpeaks)-2
		periodtrace(irespmaxpeaks(j):irespmaxpeaks(j+1)) = (irespmaxpeaks(j+1) - irespmaxpeaks(j))/fs ...
					+ (irespmaxpeaks(j+2) - 2 * irespmaxpeaks(j+1) + irespmaxpeaks(j))/fs ...
								* (((irespmaxpeaks(j):irespmaxpeaks(j+1))-irespmaxpeaks(j)) / (irespmaxpeaks(j+1) - irespmaxpeaks(j)));
	end
	periodtrace(1:irespmaxpeaks(1))=periodtrace(irespmaxpeaks(1)+1);
	periodtrace(irespmaxpeaks(end-1):end)=(irespmaxpeaks(end) - irespmaxpeaks(end-1))./fs;

	RVT = (maxtrace-mintrace)./periodtrace;
	
	RVTregressor = RVT(ivols);


	output_filename = strcat(prefix,'_RVT.txt');
	fprintf(1,'\tWriting RVT regressor to output file: %s\n',output_filename);
	fid=fopen(output_filename,'wt');
	if (RVT_td)
		RVTregressor_td = diff(RVT); RVTregressor_td = RVTregressor_td(ivols);
		if (demean); RVTregressor = RVTregressor - mean(RVTregressor); RVTregressor_td = RVTregressor_td - mean(RVTregressor_td); end
			fprintf(fid,'%f\t%f\n',[RVTregressor RVTregressor_td]');
	else
		if (demean); RVTregressor = RVTregressor - mean(RVTregressor); end
			fprintf(fid,'%f\n',RVTregressor');	
	end
	fclose(fid);

	if (RRFconv == 1)
		fprintf(1,'\tConvolving RVT regressor with RRF\n');
		t = 0:1/fs:60;
		RRF = 0.6 .* (t .^ 2.1) .* (exp(-t./1.6)) - 0.0023 .* (t .^ 3.54) .* exp(-t ./ 4.25);
		RVTconv = conv([RVT(1)*ones(60*fs,1);RVT],RRF)/sum(RRF);
		RVTconv = RVTconv(60*fs+1:60*fs+length(RVT));
		RVTconvregressor = RVTconv(ivols);
		if (RVT_td); RVTconvregressor_td = diff(RVTconv); RVTconvregressor_td = RVTconvregressor_td(ivols); end

		output_filename = strcat(prefix,'_RVT_RRFconv.txt');
		fprintf(1,'\tWriting RVT RRFconv regressor to output file: %s\n',output_filename);
		fid=fopen(output_filename,'wt');
		if (RVT_td)
			if (demean); RVTconvregressor = RVTconvregressor - mean(RVTconvregressor); RVTconvregressor_td = RVTconvregressor_td - mean(RVTconvregressor_td); end
				fprintf(fid,'%f\t%f\n',[RVTconvregressor RVTconvregressor_td]');
		else
			if (demean); RVTconvregressor = RVTconvregressor - mean(RVTconvregressor); end
				fprintf(fid,'%f\n',RVTconvregressor');
		end
		fclose(fid);
	end

	if (output_hr)
		output_filename = strcat(prefix,'_RVT_hr.txt');
		fprintf(1,'\tWriting high resolution RVT regressor to output file: %s\n',output_filename);
		fid=fopen(output_filename,'wt');
		if (RVT_td)
			if (RRFconv)
				RVTtd = diff(RVTconv); RVTtd(end+1) = RVTtd(end);
				if (demean); RVT = RVT - mean(RVT); RVTconv = RVTconv - mean(RVTconv(1:length(RVT))); RVTtd = RVTtd - mean(RVTtd(1:length(RVT))); end
					fprintf(fid,'%f\t%f\t%f\t%d\n',[RVT RVTconv RVTtd vol]');
			else
				RVTtd = diff(RVT); RVTtd(end+1) = RVTtd(end);
				if (demean); RVT = RVT - mean(RVT); RVTtd = RVTtd - mean(RVTtd); end
					fprintf(fid,'%f\t%f\t%d\n',[RVT RVTtd vol]');
			end
		else
			if (RRFconv)
				if (demean); RVT = RVT - mean(RVT); RVTconv = RVTconv - mean(RVTconv(1:length(RVT)));end
							fprintf(fid,'%f\t%f\t%d\n',[RVT RVTconv vol]');
			else
				if (demean); RVT = RVT - mean(RVT); end
					fprintf(fid,'%f\t%d\n',[RVT vol]');
			end
		end
		fclose(fid);
	end



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate and write out CO2 regressor %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (doCO2)
	fprintf(1,'Calculating CO2 regressor\n');
	
	if (exist('iCO2peaks') == 0)
		fprintf(1,'\n\nOutput of CO2 peak detection algorithm does not exist\nExiting!\n\n');
		return
	end
	
	maxenv = [];
	maxenv(1:iCO2peaks(1)) = CO2filt(iCO2peaks(1));
	for j = 1:length(iCO2peaks)-1;
		maxenv(iCO2peaks(j):iCO2peaks(j+1)) = ...
					CO2filt(iCO2peaks(j)) + (CO2filt(iCO2peaks(j+1))-CO2filt(iCO2peaks(j))) ...
								* (((iCO2peaks(j):iCO2peaks(j+1))-iCO2peaks(j)) / (iCO2peaks(j+1) - iCO2peaks(j)));
	end
	maxenv(iCO2peaks(end):length(CO2filt))=CO2filt(iCO2peaks(end));
	maxenv = maxenv';

	CO2regressor = maxenv(ivols);

	output_filename = strcat(prefix,'_CO2.txt');
	fprintf(1,'\tWriting CO2 regressor to output file: %s\n',output_filename);
	fid=fopen(output_filename,'wt');
	if (CO2_td)
		CO2regressor_td = diff(CO2); CO2regressor_td = CO2regressor_td(ivols);
		if (demean); CO2regressor = CO2regressor - mean(CO2regressor); CO2regressor_td = CO2regressor_td - mean(CO2regressor_td); end
			fprintf(fid,'%f\t%f\n',[CO2regressor CO2regressor_td]');
	else
		if (demean); CO2regressor = CO2regressor - mean(CO2regressor); end
			fprintf(fid,'%f\n',CO2regressor');	
	end
	fclose(fid);

	if (HRFconv == 1)
		fprintf(1,'\tConvolving CO2 regressor with HRF\n');
		t = 0:1/fs:25;
		HRF = exp(-t) .* ((0.00833333 .* t .^ 5) - (1.27e-13 .* t .^ 15)) ;
		CO2conv = conv([maxenv(1)*ones(60*fs,1);maxenv],HRF)/sum(HRF);
		CO2conv = CO2conv(60*fs+1:60*fs+length(maxenv));
		CO2convregressor = CO2conv(ivols);
		if (CO2_td); CO2convregressor_td = diff(CO2conv); CO2convregressor_td = CO2convregressor_td(ivols); end

		output_filename = strcat(prefix,'_CO2_HRFconv.txt');
		fprintf(1,'\tWriting CO2 HRFconv regressor to output file: %s\n',output_filename);
		fid=fopen(output_filename,'wt');
		if (CO2_td)
			if (demean); CO2convregressor = CO2convregressor - mean(CO2convregressor); CO2convregressor_td = CO2convregressor_td - mean(CO2convregressor_td); end
				fprintf(fid,'%f\t%f\n',[CO2convregressor CO2convregressor_td]');
		else
			if (demean); CO2convregressor = CO2convregressor - mean(CO2convregressor); end
				fprintf(fid,'%f\n',CO2convregressor');
		end
		fclose(fid);
	end

	if (output_hr)
		output_filename = strcat(prefix,'_CO2_hr.txt');
		fprintf(1,'\tWriting high resolution CO2 regressor to output file: %s\n',output_filename);
		fid=fopen(output_filename,'wt');
		if (CO2_td)
			if (HRFconv)
				CO2td = diff(CO2conv); CO2td(end+1) = CO2td(end);
				if (demean); maxenv = maxenv - mean(maxenv); CO2conv = CO2conv - mean(CO2conv(1:length(maxenv))); CO2td = CO2td - mean(CO2td(1:length(maxenv))); end
					fprintf(fid,'%f\t%f\t%f\t%d\n',[maxenv CO2conv CO2td vol]');
			else
				CO2td = diff(maxenv); CO2td(end+1) = CO2td(end);
				if (demean); maxenv = maxenv - mean(maxenv); CO2td = CO2td - mean(CO2td); end
					fprintf(fid,'%f\t%f\t%d\n',[maxenv CO2td vol]');
			end
		else
			if (HRFconv)
				if (demean); maxenv = maxenv - mean(maxenv); CO2conv = CO2conv - mean(CO2conv(1:length(maxenv)));end
						fprintf(fid,'%f\t%f\t%d\n',[maxenv CO2conv vol]');
			else
				if (demean); maxenv = maxenv - mean(maxenv); end
					fprintf(fid,'%f\t%d\n',[maxenv vol]');
			end
		end
		fclose(fid);
	end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate and write out O2 regressor %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (doO2)
	fprintf(1,'Calculating O2 regressor\n');

	if (exist('iO2peaks') == 0)
		fprintf(1,'\n\nOutput of O2 peak detection algorithm does not exist\nExiting!\n\n');
		return
	end


	O2env(1:iO2peaks(1))=O2filt(iO2peaks(1));
	for j=1:length(iO2peaks)-1;
		O2env(iO2peaks(j):iO2peaks(j+1))=O2filt(iO2peaks(j)) + (O2filt(iO2peaks(j+1))-O2filt(iO2peaks(j))) * (((iO2peaks(j):iO2peaks(j+1))-iO2peaks(j)) / (iO2peaks(j+1) - iO2peaks(j)));
	end
	O2env(iO2peaks(end):length(O2filt))=O2filt(iO2peaks(end));
	O2env = O2env';

	O2regressor = O2env(ivols);

	output_filename = strcat(prefix,'_O2.txt');
	fprintf(1,'\tWriting O2 regressor to output file: %s\n',output_filename);
	fid=fopen(output_filename,'wt');
	if (O2_td)
		O2regressor_td = diff(O2); O2regressor_td = O2regressor_td(ivols);
		if (demean); O2regressor = O2regressor - mean(O2regressor); O2regressor_td = O2regressor_td - mean(O2regressor_td); end
			fprintf(fid,'%f\t%f\n',[O2regressor O2regressor_td]');
	else
		if (demean); O2regressor = O2regressor - mean(O2regressor); end
			fprintf(fid,'%f\n',O2regressor');	
	end
	fclose(fid);

	if (HRFconv == 1)
		fprintf(1,'\tConvolving O2 regressor with HRF\n');
		t = 0:1/fs:25;
		HRF = exp(-t) .* ((0.00833333 .* t .^ 5) - (1.27e-13 .* t .^ 15)) ;
		O2conv = conv([O2env(1)*ones(60*fs,1);O2env],HRF)/sum(HRF);
		O2conv = O2conv(60*fs+1:60*fs+length(O2env));
		O2convregressor = O2conv(ivols);
		if (O2_td); O2convregressor_td = diff(O2conv); O2convregressor_td = O2convregressor_td(ivols); end

		output_filename = strcat(prefix,'_O2_HRFconv.txt');
		fprintf(1,'\tWriting O2 HRFconv regressor to output file: %s\n',output_filename);
		fid=fopen(output_filename,'wt');
		if (O2_td)
			if (demean); O2convregressor = O2convregressor - mean(O2convregressor); O2convregressor_td = O2convregressor_td - mean(O2convregressor_td); end
				fprintf(fid,'%f\t%f\n',[O2convregressor O2convregressor_td]');
		else
			if (demean); O2convregressor = O2convregressor - mean(O2convregressor); end
				fprintf(fid,'%f\n',O2convregressor');
		end
		fclose(fid);
	end

	if (output_hr)
		output_filename = strcat(prefix,'_O2_hr.txt');
		fprintf(1,'\tWriting high resolution O2 regressor to output file: %s\n',output_filename);
		fid=fopen(output_filename,'wt');
		if (O2_td)
			if (HRFconv)
				O2td = diff(O2conv); O2td(end+1) = O2td(end);
				if (demean); O2env = O2env - mean(O2env); O2conv = O2conv - mean(O2conv(1:length(O2env))); O2td = O2td - mean(O2td(1:length(O2env))); end
					fprintf(fid,'%f\t%f\t%f\t%d\n',[O2env O2conv O2td vol]');
			else
				O2td = diff(O2env); O2td(end+1) = O2td(end);
				if (demean); O2env = O2env - mean(O2env); O2td = O2td - mean(O2td); end
					fprintf(fid,'%f\t%f\t%d\n',[O2env O2td vol]');
			end
		else
			if (HRFconv)
				if (demean); O2env = O2env - mean(O2env); O2conv = O2conv - mean(O2conv(1:length(O2env)));end
						fprintf(fid,'%f\t%f\t%d\n',[O2env O2conv vol]');
			else
				if (demean); O2env = O2env - mean(O2env); end
					fprintf(fid,'%f\t%d\n',[O2env vol]');
			end
		end
		fclose(fid);
	end
end


fprintf(1,'\nFinished!\n\n');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% End of calc_phys_regressor script %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% The following subfunctions are included so that this file will  %%%%%%%%
%%%%%%%% work as a standalone matlab script. This functions have been    %%%%%%%%
%%%%%%%% borrowed from Jimmy Shen (jimmy@rotman-baycrest.on.ca). See     %%%%%%%%
%%%%%%%% http://www.rotman-baycrest.on.ca/~jimmy/NIFTI/ for more details %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2009, Jimmy Shen                                                %
% All rights reserved.                                                          %
%                                                                               %
% Redistribution and use in source and binary forms, with or without            % 
% modification, are permitted provided that the following conditions are        % 
% met:                                                                          % 
%                                                                               % 
%     * Redistributions of source code must retain the above copyright          % 
%       notice, this list of conditions and the following disclaimer.           % 
%     * Redistributions in binary form must reproduce the above copyright       % 
%       notice, this list of conditions and the following disclaimer in         % 
%       the documentation and/or other materials provided with the distribution % 
%                                                                               % 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"   %
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE     %                                                                               % 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE    % 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE      % 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR           % 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF          % 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS      % 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN       %                                                                              % 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)       % 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE    % 
% POSSIBILITY OF SUCH DAMAGE.                                                   % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function nii = load_untouch_nii(filename, img_idx, dim5_idx, dim6_idx, dim7_idx, ...
			old_RGB, slice_idx)

   if ~exist('filename','var')
      error('Usage: nii = load_untouch_nii(filename, [img_idx], [dim5_idx], [dim6_idx], [dim7_idx], [old_RGB], [slice_idx])');
   end

   if ~exist('img_idx','var') | isempty(img_idx)
      img_idx = [];
   end

   if ~exist('dim5_idx','var') | isempty(dim5_idx)
      dim5_idx = [];
   end

   if ~exist('dim6_idx','var') | isempty(dim6_idx)
      dim6_idx = [];
   end

   if ~exist('dim7_idx','var') | isempty(dim7_idx)
      dim7_idx = [];
   end

   if ~exist('old_RGB','var') | isempty(old_RGB)
      old_RGB = 0;
   end

   if ~exist('slice_idx','var') | isempty(slice_idx)
      slice_idx = [];
   end

   %  Read the dataset header
   %
   [nii.hdr,nii.filetype,nii.fileprefix,nii.machine] = load_nii_hdr(filename);

   if nii.filetype == 0
      nii.hdr = load_untouch0_nii_hdr(nii.fileprefix,nii.machine);
      nii.ext = [];
   else
      nii.hdr = load_untouch_nii_hdr(nii.fileprefix,nii.machine,nii.filetype);

      %  Read the header extension
      %
      nii.ext = load_nii_ext(filename);
   end

   %  Read the dataset body
   %
   [nii.img,nii.hdr] = load_untouch_nii_img(nii.hdr,nii.filetype,nii.fileprefix, ...
		nii.machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB,slice_idx);

   %  Perform some of sform/qform transform
   %
%   nii = xform_nii(nii, tolerance, preferredForm);

   nii.untouch = 1;

   return					% load_untouch_nii

%---------------------------------------------------------------------

function [hdr, filetype, fileprefix, machine] = load_nii_hdr(fileprefix)

   if ~exist('fileprefix','var'),
      error('Usage: [hdr, filetype, fileprefix, machine] = load_nii_hdr(filename)');
   end

   machine = 'ieee-le';
   new_ext = 0;

   if findstr('.nii',fileprefix)
      new_ext = 1;
      fileprefix = strrep(fileprefix,'.nii','');
   end

   if findstr('.hdr',fileprefix)
      fileprefix = strrep(fileprefix,'.hdr','');
   end

   if findstr('.img',fileprefix)
      fileprefix = strrep(fileprefix,'.img','');
   end

   if new_ext
      fn = sprintf('%s.nii',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.nii".', fileprefix);
         error(msg);
      end
   else
      fn = sprintf('%s.hdr',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.hdr".', fileprefix);
         error(msg);
      end
   end

   fid = fopen(fn,'r',machine);
    
   if fid < 0,
      msg = sprintf('Cannot open file %s.',fn);
      error(msg);
   else
      fseek(fid,0,'bof');

      if fread(fid,1,'int32') == 348
         hdr = read_header(fid);
         fclose(fid);
      else
         fclose(fid);

         %  first try reading the opposite endian to 'machine'
         %
         switch machine,
         case 'ieee-le', machine = 'ieee-be';
         case 'ieee-be', machine = 'ieee-le';
         end

         fid = fopen(fn,'r',machine);

         if fid < 0,
            msg = sprintf('Cannot open file %s.',fn);
            error(msg);
         else
            fseek(fid,0,'bof');

            if fread(fid,1,'int32') ~= 348

               %  Now throw an error
               %
               msg = sprintf('File "%s" is corrupted.',fn);
               error(msg);
            end

            hdr = read_header(fid);
            fclose(fid);
         end
      end
   end

   if strcmp(hdr.hist.magic, 'n+1')
      filetype = 2;
   elseif strcmp(hdr.hist.magic, 'ni1')
      filetype = 1;
   else
      filetype = 0;
   end

   return					% load_nii_hdr


%---------------------------------------------------------------------
function [ dsr ] = read_header(fid)

        %  Original header structures
	%  struct dsr
	%       { 
	%       struct header_key hk;            /*   0 +  40       */
	%       struct image_dimension dime;     /*  40 + 108       */
	%       struct data_history hist;        /* 148 + 200       */
	%       };                               /* total= 348 bytes*/

    dsr.hk   = header_key(fid);
    dsr.dime = image_dimension(fid);
    dsr.hist = data_history(fid);

    %  For Analyze data format
    %
    if ~strcmp(dsr.hist.magic, 'n+1') & ~strcmp(dsr.hist.magic, 'ni1')
        dsr.hist.qform_code = 0;
        dsr.hist.sform_code = 0;
    end

    return					% read_header


%---------------------------------------------------------------------
function [ hk ] = header_key(fid)

    fseek(fid,0,'bof');
    
	%  Original header structures	
	%  struct header_key                     /* header key      */ 
	%       {                                /* off + size      */
	%       int sizeof_hdr                   /*  0 +  4         */
	%       char data_type[10];              /*  4 + 10         */
	%       char db_name[18];                /* 14 + 18         */
	%       int extents;                     /* 32 +  4         */
	%       short int session_error;         /* 36 +  2         */
	%       char regular;                    /* 38 +  1         */
	%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
	%       };                               /* total=40 bytes  */
	%
	% int sizeof_header   Should be 348.
	% char regular        Must be 'r' to indicate that all images and 
	%                     volumes are the same size. 

    v6 = version;
    if str2num(v6(1))<6
       directchar = '*char';
    else
       directchar = 'uchar=>char';
    end

    hk.sizeof_hdr    = fread(fid, 1,'int32')';	% should be 348!
    hk.data_type     = deblank(fread(fid,10,directchar)');
    hk.db_name       = deblank(fread(fid,18,directchar)');
    hk.extents       = fread(fid, 1,'int32')';
    hk.session_error = fread(fid, 1,'int16')';
    hk.regular       = fread(fid, 1,directchar)';
    hk.dim_info      = fread(fid, 1,'uchar')';
    
    return					% header_key


%---------------------------------------------------------------------
function [ dime ] = image_dimension(fid)

	%  Original header structures    
	%  struct image_dimension
	%       {                                /* off + size      */
	%       short int dim[8];                /* 0 + 16          */
        %       /*
        %           dim[0]      Number of dimensions in database; usually 4. 
        %           dim[1]      Image X dimension;  number of *pixels* in an image row. 
        %           dim[2]      Image Y dimension;  number of *pixel rows* in slice. 
        %           dim[3]      Volume Z dimension; number of *slices* in a volume. 
        %           dim[4]      Time points; number of volumes in database
        %       */
	%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
	%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
	%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
	%       short int intent_code;   % short int unused1;   /* 28 + 2 */
	%       short int datatype;              /* 30 + 2          */
	%       short int bitpix;                /* 32 + 2          */
	%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
	%       float pixdim[8];                 /* 36 + 32         */
	%	/*
	%		pixdim[] specifies the voxel dimensions:
	%		pixdim[1] - voxel width, mm
	%		pixdim[2] - voxel height, mm
	%		pixdim[3] - slice thickness, mm
	%		pixdim[4] - volume timing, in msec
	%					..etc
	%	*/
	%       float vox_offset;                /* 68 + 4          */
	%       float scl_slope;   % float roi_scale;     /* 72 + 4 */
	%       float scl_inter;   % float funused1;      /* 76 + 4 */
	%       short slice_end;   % float funused2;      /* 80 + 2 */
	%       char slice_code;   % float funused2;      /* 82 + 1 */
	%       char xyzt_units;   % float funused2;      /* 83 + 1 */
	%       float cal_max;                   /* 84 + 4          */
	%       float cal_min;                   /* 88 + 4          */
	%       float slice_duration;   % int compressed; /* 92 + 4 */
	%       float toffset;   % int verified;          /* 96 + 4 */
	%       int glmax;                       /* 100 + 4         */
	%       int glmin;                       /* 104 + 4         */
	%       };                               /* total=108 bytes */
	
    dime.dim        = fread(fid,8,'int16')';
    dime.intent_p1  = fread(fid,1,'float32')';
    dime.intent_p2  = fread(fid,1,'float32')';
    dime.intent_p3  = fread(fid,1,'float32')';
    dime.intent_code = fread(fid,1,'int16')';
    dime.datatype   = fread(fid,1,'int16')';
    dime.bitpix     = fread(fid,1,'int16')';
    dime.slice_start = fread(fid,1,'int16')';
    dime.pixdim     = fread(fid,8,'float32')';
    dime.vox_offset = fread(fid,1,'float32')';
    dime.scl_slope  = fread(fid,1,'float32')';
    dime.scl_inter  = fread(fid,1,'float32')';
    dime.slice_end  = fread(fid,1,'int16')';
    dime.slice_code = fread(fid,1,'uchar')';
    dime.xyzt_units = fread(fid,1,'uchar')';
    dime.cal_max    = fread(fid,1,'float32')';
    dime.cal_min    = fread(fid,1,'float32')';
    dime.slice_duration = fread(fid,1,'float32')';
    dime.toffset    = fread(fid,1,'float32')';
    dime.glmax      = fread(fid,1,'int32')';
    dime.glmin      = fread(fid,1,'int32')';
        
    return					% image_dimension


%---------------------------------------------------------------------
function [ hist ] = data_history(fid)
        
	%  Original header structures
	%  struct data_history       
	%       {                                /* off + size      */
	%       char descrip[80];                /* 0 + 80          */
	%       char aux_file[24];               /* 80 + 24         */
	%       short int qform_code;            /* 104 + 2         */
	%       short int sform_code;            /* 106 + 2         */
	%       float quatern_b;                 /* 108 + 4         */
	%       float quatern_c;                 /* 112 + 4         */
	%       float quatern_d;                 /* 116 + 4         */
	%       float qoffset_x;                 /* 120 + 4         */
	%       float qoffset_y;                 /* 124 + 4         */
	%       float qoffset_z;                 /* 128 + 4         */
	%       float srow_x[4];                 /* 132 + 16        */
	%       float srow_y[4];                 /* 148 + 16        */
	%       float srow_z[4];                 /* 164 + 16        */
	%       char intent_name[16];            /* 180 + 16        */
	%       char magic[4];   % int smin;     /* 196 + 4         */
	%       };                               /* total=200 bytes */

    v6 = version;
    if str2num(v6(1))<6
       directchar = '*char';
    else
       directchar = 'uchar=>char';
    end

    hist.descrip     = deblank(fread(fid,80,directchar)');
    hist.aux_file    = deblank(fread(fid,24,directchar)');
    hist.qform_code  = fread(fid,1,'int16')';
    hist.sform_code  = fread(fid,1,'int16')';
    hist.quatern_b   = fread(fid,1,'float32')';
    hist.quatern_c   = fread(fid,1,'float32')';
    hist.quatern_d   = fread(fid,1,'float32')';
    hist.qoffset_x   = fread(fid,1,'float32')';
    hist.qoffset_y   = fread(fid,1,'float32')';
    hist.qoffset_z   = fread(fid,1,'float32')';
    hist.srow_x      = fread(fid,4,'float32')';
    hist.srow_y      = fread(fid,4,'float32')';
    hist.srow_z      = fread(fid,4,'float32')';
    hist.intent_name = deblank(fread(fid,16,directchar)');
    hist.magic       = deblank(fread(fid,4,directchar)');

    fseek(fid,253,'bof');
    hist.originator  = fread(fid, 5,'int16')';
    
    return					% data_history


%---------------------------------------------------------------------

function hdr = load_untouch_nii_hdr(fileprefix, machine, filetype)

   if filetype == 2
      fn = sprintf('%s.nii',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.nii".', fileprefix);
         error(msg);
      end
   else
      fn = sprintf('%s.hdr',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.hdr".', fileprefix);
         error(msg);
      end
   end

   fid = fopen(fn,'r',machine);
    
   if fid < 0,
      msg = sprintf('Cannot open file %s.',fn);
      error(msg);
   else
      fseek(fid,0,'bof');
      hdr = read_header1(fid);
      fclose(fid);
   end

   return					% load_nii_hdr


%---------------------------------------------------------------------
function [ dsr ] = read_header1(fid)

        %  Original header structures
	%  struct dsr
	%       { 
	%       struct header_key hk;            /*   0 +  40       */
	%       struct image_dimension dime;     /*  40 + 108       */
	%       struct data_history hist;        /* 148 + 200       */
	%       };                               /* total= 348 bytes*/

    dsr.hk   = header_key1(fid);
    dsr.dime = image_dimension1(fid);
    dsr.hist = data_history1(fid);

    %  For Analyze data format
    %
    if ~strcmp(dsr.hist.magic, 'n+1') & ~strcmp(dsr.hist.magic, 'ni1')
        dsr.hist.qform_code = 0;
        dsr.hist.sform_code = 0;
    end

    return					% read_header


%---------------------------------------------------------------------
function [ hk ] = header_key1(fid)

    fseek(fid,0,'bof');
    
	%  Original header structures	
	%  struct header_key                     /* header key      */ 
	%       {                                /* off + size      */
	%       int sizeof_hdr                   /*  0 +  4         */
	%       char data_type[10];              /*  4 + 10         */
	%       char db_name[18];                /* 14 + 18         */
	%       int extents;                     /* 32 +  4         */
	%       short int session_error;         /* 36 +  2         */
	%       char regular;                    /* 38 +  1         */
	%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
	%       };                               /* total=40 bytes  */
	%
	% int sizeof_header   Should be 348.
	% char regular        Must be 'r' to indicate that all images and 
	%                     volumes are the same size. 

    v6 = version;
    if str2num(v6(1))<6
       directchar = '*char';
    else
       directchar = 'uchar=>char';
    end
	
    hk.sizeof_hdr    = fread(fid, 1,'int32')';	% should be 348!
    hk.data_type     = deblank(fread(fid,10,directchar)');
    hk.db_name       = deblank(fread(fid,18,directchar)');
    hk.extents       = fread(fid, 1,'int32')';
    hk.session_error = fread(fid, 1,'int16')';
    hk.regular       = fread(fid, 1,directchar)';
    hk.dim_info      = fread(fid, 1,'uchar')';
    
    return					% header_key


%---------------------------------------------------------------------
function [ dime ] = image_dimension1(fid)

	%  Original header structures    
	%  struct image_dimension
	%       {                                /* off + size      */
	%       short int dim[8];                /* 0 + 16          */
        %       /*
        %           dim[0]      Number of dimensions in database; usually 4. 
        %           dim[1]      Image X dimension;  number of *pixels* in an image row. 
        %           dim[2]      Image Y dimension;  number of *pixel rows* in slice. 
        %           dim[3]      Volume Z dimension; number of *slices* in a volume. 
        %           dim[4]      Time points; number of volumes in database
        %       */
	%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
	%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
	%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
	%       short int intent_code;   % short int unused1;   /* 28 + 2 */
	%       short int datatype;              /* 30 + 2          */
	%       short int bitpix;                /* 32 + 2          */
	%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
	%       float pixdim[8];                 /* 36 + 32         */
	%	/*
	%		pixdim[] specifies the voxel dimensions:
	%		pixdim[1] - voxel width, mm
	%		pixdim[2] - voxel height, mm
	%		pixdim[3] - slice thickness, mm
	%		pixdim[4] - volume timing, in msec
	%					..etc
	%	*/
	%       float vox_offset;                /* 68 + 4          */
	%       float scl_slope;   % float roi_scale;     /* 72 + 4 */
	%       float scl_inter;   % float funused1;      /* 76 + 4 */
	%       short slice_end;   % float funused2;      /* 80 + 2 */
	%       char slice_code;   % float funused2;      /* 82 + 1 */
	%       char xyzt_units;   % float funused2;      /* 83 + 1 */
	%       float cal_max;                   /* 84 + 4          */
	%       float cal_min;                   /* 88 + 4          */
	%       float slice_duration;   % int compressed; /* 92 + 4 */
	%       float toffset;   % int verified;          /* 96 + 4 */
	%       int glmax;                       /* 100 + 4         */
	%       int glmin;                       /* 104 + 4         */
	%       };                               /* total=108 bytes */
	
    dime.dim        = fread(fid,8,'int16')';
    dime.intent_p1  = fread(fid,1,'float32')';
    dime.intent_p2  = fread(fid,1,'float32')';
    dime.intent_p3  = fread(fid,1,'float32')';
    dime.intent_code = fread(fid,1,'int16')';
    dime.datatype   = fread(fid,1,'int16')';
    dime.bitpix     = fread(fid,1,'int16')';
    dime.slice_start = fread(fid,1,'int16')';
    dime.pixdim     = fread(fid,8,'float32')';
    dime.vox_offset = fread(fid,1,'float32')';
    dime.scl_slope  = fread(fid,1,'float32')';
    dime.scl_inter  = fread(fid,1,'float32')';
    dime.slice_end  = fread(fid,1,'int16')';
    dime.slice_code = fread(fid,1,'uchar')';
    dime.xyzt_units = fread(fid,1,'uchar')';
    dime.cal_max    = fread(fid,1,'float32')';
    dime.cal_min    = fread(fid,1,'float32')';
    dime.slice_duration = fread(fid,1,'float32')';
    dime.toffset    = fread(fid,1,'float32')';
    dime.glmax      = fread(fid,1,'int32')';
    dime.glmin      = fread(fid,1,'int32')';
        
    return					% image_dimension


%---------------------------------------------------------------------
function [ hist ] = data_history1(fid)
        
	%  Original header structures
	%  struct data_history       
	%       {                                /* off + size      */
	%       char descrip[80];                /* 0 + 80          */
	%       char aux_file[24];               /* 80 + 24         */
	%       short int qform_code;            /* 104 + 2         */
	%       short int sform_code;            /* 106 + 2         */
	%       float quatern_b;                 /* 108 + 4         */
	%       float quatern_c;                 /* 112 + 4         */
	%       float quatern_d;                 /* 116 + 4         */
	%       float qoffset_x;                 /* 120 + 4         */
	%       float qoffset_y;                 /* 124 + 4         */
	%       float qoffset_z;                 /* 128 + 4         */
	%       float srow_x[4];                 /* 132 + 16        */
	%       float srow_y[4];                 /* 148 + 16        */
	%       float srow_z[4];                 /* 164 + 16        */
	%       char intent_name[16];            /* 180 + 16        */
	%       char magic[4];   % int smin;     /* 196 + 4         */
	%       };                               /* total=200 bytes */

    v6 = version;
    if str2num(v6(1))<6
       directchar = '*char';
    else
       directchar = 'uchar=>char';
    end
    
    hist.descrip     = deblank(fread(fid,80,directchar)');
    hist.aux_file    = deblank(fread(fid,24,directchar)');
    hist.qform_code  = fread(fid,1,'int16')';
    hist.sform_code  = fread(fid,1,'int16')';
    hist.quatern_b   = fread(fid,1,'float32')';
    hist.quatern_c   = fread(fid,1,'float32')';
    hist.quatern_d   = fread(fid,1,'float32')';
    hist.qoffset_x   = fread(fid,1,'float32')';
    hist.qoffset_y   = fread(fid,1,'float32')';
    hist.qoffset_z   = fread(fid,1,'float32')';
    hist.srow_x      = fread(fid,4,'float32')';
    hist.srow_y      = fread(fid,4,'float32')';
    hist.srow_z      = fread(fid,4,'float32')';
    hist.intent_name = deblank(fread(fid,16,directchar)');
    hist.magic       = deblank(fread(fid,4,directchar)');
    
    return					% data_history

%---------------------------------------------------------------------

function ext = load_nii_ext(fileprefix)

   if ~exist('fileprefix','var'),
      error('Usage: ext = load_nii_ext(filename)');
   end

   machine = 'ieee-le';
   new_ext = 0;

   if findstr('.nii',fileprefix)
      new_ext = 1;
      fileprefix = strrep(fileprefix,'.nii','');
   end

   if findstr('.hdr',fileprefix)
      fileprefix = strrep(fileprefix,'.hdr','');
   end

   if findstr('.img',fileprefix)
      fileprefix = strrep(fileprefix,'.img','');
   end

   if new_ext
      fn = sprintf('%s.nii',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.nii".', fileprefix);
         error(msg);
      end
   else
      fn = sprintf('%s.hdr',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.hdr".', fileprefix);
         error(msg);
      end
   end

   fid = fopen(fn,'r',machine);
   vox_offset = 0;
    
   if fid < 0,
      msg = sprintf('Cannot open file %s.',fn);
      error(msg);
   else
      fseek(fid,0,'bof');

      if fread(fid,1,'int32') == 348
         if new_ext
            fseek(fid,108,'bof');
            vox_offset = fread(fid,1,'float32');
         end

         ext = read_extension(fid, vox_offset);
         fclose(fid);
      else
         fclose(fid);

         %  first try reading the opposite endian to 'machine'
         %
         switch machine,
         case 'ieee-le', machine = 'ieee-be';
         case 'ieee-be', machine = 'ieee-le';
         end

         fid = fopen(fn,'r',machine);

         if fid < 0,
            msg = sprintf('Cannot open file %s.',fn);
            error(msg);
         else
            fseek(fid,0,'bof');

            if fread(fid,1,'int32') ~= 348

               %  Now throw an error
               %
               msg = sprintf('File "%s" is corrupted.',fn);
               error(msg);
            end

            if new_ext
               fseek(fid,108,'bof');
               vox_offset = fread(fid,1,'float32');
            end

            ext = read_extension(fid, vox_offset);
            fclose(fid);
         end
      end
   end

   return                                       % load_nii_ext


%---------------------------------------------------------------------
function ext = read_extension(fid, vox_offset)

   ext = [];

   if vox_offset
      end_of_ext = vox_offset;
   else
      fseek(fid, 0, 'eof');
      end_of_ext = ftell(fid);
   end

   if end_of_ext > 352
      fseek(fid, 348, 'bof');
      ext.extension = fread(fid,4)';
   end

   if isempty(ext) | ext.extension(1) == 0
      ext = [];
      return;
   end

   i = 1;

   while(ftell(fid) < end_of_ext)
      ext.section(i).esize = fread(fid,1,'int32');
      ext.section(i).ecode = fread(fid,1,'int32');
      ext.section(i).edata = char(fread(fid,ext.section(i).esize-8)');
      i = i + 1;
   end

   ext.num_ext = length(ext.section);

   return                                               % read_extension
%---------------------------------------------------------------------

function [img,hdr] = load_untouch_nii_img(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB,slice_idx)

   if ~exist('hdr','var') | ~exist('filetype','var') | ~exist('fileprefix','var') | ~exist('machine','var')
      error('Usage: [img,hdr] = load_nii_img(hdr,filetype,fileprefix,machine,[img_idx],[dim5_idx],[dim6_idx],[dim7_idx],[old_RGB],[slice_idx]);');
   end

   if ~exist('img_idx','var') | isempty(img_idx) | hdr.dime.dim(5)<1
      img_idx = [];
   end

   if ~exist('dim5_idx','var') | isempty(dim5_idx) | hdr.dime.dim(6)<1
      dim5_idx = [];
   end

   if ~exist('dim6_idx','var') | isempty(dim6_idx) | hdr.dime.dim(7)<1
      dim6_idx = [];
   end

   if ~exist('dim7_idx','var') | isempty(dim7_idx) | hdr.dime.dim(8)<1
      dim7_idx = [];
   end

   if ~exist('old_RGB','var') | isempty(old_RGB)
      old_RGB = 0;
   end

   if ~exist('slice_idx','var') | isempty(slice_idx) | hdr.dime.dim(4)<1
      slice_idx = [];
   end

   %  check img_idx
   %
   if ~isempty(img_idx) & ~isnumeric(img_idx)
      error('"img_idx" should be a numerical array.');
   end

   if length(unique(img_idx)) ~= length(img_idx)
      error('Duplicate image index in "img_idx"');
   end

   if ~isempty(img_idx) & (min(img_idx) < 1 | max(img_idx) > hdr.dime.dim(5))
      max_range = hdr.dime.dim(5);

      if max_range == 1
         error(['"img_idx" should be 1.']);
      else
         range = ['1 ' num2str(max_range)];
         error(['"img_idx" should be an integer within the range of [' range '].']);
      end
   end

   %  check dim5_idx
   %
   if ~isempty(dim5_idx) & ~isnumeric(dim5_idx)
      error('"dim5_idx" should be a numerical array.');
   end

   if length(unique(dim5_idx)) ~= length(dim5_idx)
      error('Duplicate index in "dim5_idx"');
   end

   if ~isempty(dim5_idx) & (min(dim5_idx) < 1 | max(dim5_idx) > hdr.dime.dim(6))
      max_range = hdr.dime.dim(6);

      if max_range == 1
         error(['"dim5_idx" should be 1.']);
      else
         range = ['1 ' num2str(max_range)];
         error(['"dim5_idx" should be an integer within the range of [' range '].']);
      end
   end

   %  check dim6_idx
   %
   if ~isempty(dim6_idx) & ~isnumeric(dim6_idx)
      error('"dim6_idx" should be a numerical array.');
   end

   if length(unique(dim6_idx)) ~= length(dim6_idx)
      error('Duplicate index in "dim6_idx"');
   end

   if ~isempty(dim6_idx) & (min(dim6_idx) < 1 | max(dim6_idx) > hdr.dime.dim(7))
      max_range = hdr.dime.dim(7);

      if max_range == 1
         error(['"dim6_idx" should be 1.']);
      else
         range = ['1 ' num2str(max_range)];
         error(['"dim6_idx" should be an integer within the range of [' range '].']);
      end
   end

   %  check dim7_idx
   %
   if ~isempty(dim7_idx) & ~isnumeric(dim7_idx)
      error('"dim7_idx" should be a numerical array.');
   end

   if length(unique(dim7_idx)) ~= length(dim7_idx)
      error('Duplicate index in "dim7_idx"');
   end

   if ~isempty(dim7_idx) & (min(dim7_idx) < 1 | max(dim7_idx) > hdr.dime.dim(8))
      max_range = hdr.dime.dim(8);

      if max_range == 1
         error(['"dim7_idx" should be 1.']);
      else
         range = ['1 ' num2str(max_range)];
         error(['"dim7_idx" should be an integer within the range of [' range '].']);
      end
   end

   %  check slice_idx
   %
   if ~isempty(slice_idx) & ~isnumeric(slice_idx)
      error('"slice_idx" should be a numerical array.');
   end

   if length(unique(slice_idx)) ~= length(slice_idx)
      error('Duplicate index in "slice_idx"');
   end

   if ~isempty(slice_idx) & (min(slice_idx) < 1 | max(slice_idx) > hdr.dime.dim(4))
      max_range = hdr.dime.dim(4);

      if max_range == 1
         error(['"slice_idx" should be 1.']);
      else
         range = ['1 ' num2str(max_range)];
         error(['"slice_idx" should be an integer within the range of [' range '].']);
      end
   end

   [img,hdr] = read_image(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB,slice_idx);

   return					% load_nii_img


%---------------------------------------------------------------------
function [img,hdr] = read_image(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB,slice_idx)

   switch filetype
   case {0, 1}
      fn = [fileprefix '.img'];
   case 2
      fn = [fileprefix '.nii'];
   end

   fid = fopen(fn,'r',machine);

   if fid < 0,
      msg = sprintf('Cannot open file %s.',fn);
      error(msg);
   end

   %  Set bitpix according to datatype
   %
   %  /*Acceptable values for datatype are*/ 
   %
   %     0 None                     (Unknown bit per voxel) % DT_NONE, DT_UNKNOWN 
   %     1 Binary                         (ubit1, bitpix=1) % DT_BINARY 
   %     2 Unsigned char         (uchar or uint8, bitpix=8) % DT_UINT8, NIFTI_TYPE_UINT8 
   %     4 Signed short                  (int16, bitpix=16) % DT_INT16, NIFTI_TYPE_INT16 
   %     8 Signed integer                (int32, bitpix=32) % DT_INT32, NIFTI_TYPE_INT32 
   %    16 Floating point    (single or float32, bitpix=32) % DT_FLOAT32, NIFTI_TYPE_FLOAT32 
   %    32 Complex, 2 float32      (Use float32, bitpix=64) % DT_COMPLEX64, NIFTI_TYPE_COMPLEX64
   %    64 Double precision  (double or float64, bitpix=64) % DT_FLOAT64, NIFTI_TYPE_FLOAT64 
   %   128 uint8 RGB                 (Use uint8, bitpix=24) % DT_RGB24, NIFTI_TYPE_RGB24 
   %   256 Signed char            (schar or int8, bitpix=8) % DT_INT8, NIFTI_TYPE_INT8 
   %   511 Single RGB              (Use float32, bitpix=96) % DT_RGB96, NIFTI_TYPE_RGB96
   %   512 Unsigned short               (uint16, bitpix=16) % DT_UNINT16, NIFTI_TYPE_UNINT16 
   %   768 Unsigned integer             (uint32, bitpix=32) % DT_UNINT32, NIFTI_TYPE_UNINT32 
   %  1024 Signed long long              (int64, bitpix=64) % DT_INT64, NIFTI_TYPE_INT64
   %  1280 Unsigned long long           (uint64, bitpix=64) % DT_UINT64, NIFTI_TYPE_UINT64 
   %  1536 Long double, float128  (Unsupported, bitpix=128) % DT_FLOAT128, NIFTI_TYPE_FLOAT128 
   %  1792 Complex128, 2 float64  (Use float64, bitpix=128) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
   %  2048 Complex256, 2 float128 (Unsupported, bitpix=256) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
   %
   switch hdr.dime.datatype
   case   1,
      hdr.dime.bitpix = 1;  precision = 'ubit1';
   case   2,
      hdr.dime.bitpix = 8;  precision = 'uint8';
   case   4,
      hdr.dime.bitpix = 16; precision = 'int16';
   case   8,
      hdr.dime.bitpix = 32; precision = 'int32';
   case  16,
      hdr.dime.bitpix = 32; precision = 'float32';
   case  32,
      hdr.dime.bitpix = 64; precision = 'float32';
   case  64,
      hdr.dime.bitpix = 64; precision = 'float64';
   case 128,
      hdr.dime.bitpix = 24; precision = 'uint8';
   case 256 
      hdr.dime.bitpix = 8;  precision = 'int8';
   case 511 
      hdr.dime.bitpix = 96; precision = 'float32';
   case 512 
      hdr.dime.bitpix = 16; precision = 'uint16';
   case 768 
      hdr.dime.bitpix = 32; precision = 'uint32';
   case 1024
      hdr.dime.bitpix = 64; precision = 'int64';
   case 1280
      hdr.dime.bitpix = 64; precision = 'uint64';
   case 1792,
      hdr.dime.bitpix = 128; precision = 'float64';
   otherwise
      error('This datatype is not supported'); 
   end

   tmp = hdr.dime.dim(2:end);
   tmp(find(tmp < 1)) = 1;
   hdr.dime.dim(2:end) = tmp;

   %  move pointer to the start of image block
   %
   switch filetype
   case {0, 1}
      fseek(fid, 0, 'bof');
   case 2
      fseek(fid, hdr.dime.vox_offset, 'bof');
   end

   %  Load whole image block for old Analyze format or binary image;
   %  otherwise, load images that are specified in img_idx, dim5_idx,
   %  dim6_idx, and dim7_idx
   %
   %  For binary image, we have to read all because pos can not be
   %  seeked in bit and can not be calculated the way below.
   %
   if hdr.dime.datatype == 1 | isequal(hdr.dime.dim(4:8),ones(1,5)) | ...
	(isempty(img_idx) & isempty(dim5_idx) & isempty(dim6_idx) & isempty(dim7_idx) & isempty(slice_idx))

      %  For each frame, precision of value will be read 
      %  in img_siz times, where img_siz is only the 
      %  dimension size of an image, not the byte storage
      %  size of an image.
      %
      img_siz = prod(hdr.dime.dim(2:8));

      %  For complex float32 or complex float64, voxel values
      %  include [real, imag]
      %
      if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
         img_siz = img_siz * 2;
      end
	 
      %MPH: For RGB24, voxel values include 3 separate color planes
      %
      if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
	 img_siz = img_siz * 3;
      end

      img = fread(fid, img_siz, sprintf('*%s',precision));

      d1 = hdr.dime.dim(2);
      d2 = hdr.dime.dim(3);
      d3 = hdr.dime.dim(4);
      d4 = hdr.dime.dim(5);
      d5 = hdr.dime.dim(6);
      d6 = hdr.dime.dim(7);
      d7 = hdr.dime.dim(8);

      if isempty(slice_idx)
         slice_idx = 1:d3;
      end

      if isempty(img_idx)
         img_idx = 1:d4;
      end

      if isempty(dim5_idx)
         dim5_idx = 1:d5;
      end

      if isempty(dim6_idx)
         dim6_idx = 1:d6;
      end

      if isempty(dim7_idx)
         dim7_idx = 1:d7;
      end
   else

      d1 = hdr.dime.dim(2);
      d2 = hdr.dime.dim(3);
      d3 = hdr.dime.dim(4);
      d4 = hdr.dime.dim(5);
      d5 = hdr.dime.dim(6);
      d6 = hdr.dime.dim(7);
      d7 = hdr.dime.dim(8);

      if isempty(slice_idx)
         slice_idx = 1:d3;
      end

      if isempty(img_idx)
         img_idx = 1:d4;
      end

      if isempty(dim5_idx)
         dim5_idx = 1:d5;
      end

      if isempty(dim6_idx)
         dim6_idx = 1:d6;
      end

      if isempty(dim7_idx)
         dim7_idx = 1:d7;
      end
      
      %ROMAN: begin
      roman = 1;
      if(roman)

         %  compute size of one slice
         %
         img_siz = prod(hdr.dime.dim(2:3));

         %  For complex float32 or complex float64, voxel values
         %  include [real, imag]
         %
         if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
            img_siz = img_siz * 2;
         end

         %MPH: For RGB24, voxel values include 3 separate color planes
         %
         if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
            img_siz = img_siz * 3;
         end

         % preallocate img
         img = zeros(img_siz, length(slice_idx)*length(img_idx)*length(dim5_idx)*length(dim6_idx)*length(dim7_idx) );
         currentIndex = 1;
      else
        img = [];
      end; %if(roman)
      % ROMAN: end

      for i7=1:length(dim7_idx)
         for i6=1:length(dim6_idx)
            for i5=1:length(dim5_idx)
               for t=1:length(img_idx)
               for s=1:length(slice_idx)

                  %  Position is seeked in bytes. To convert dimension size
                  %  to byte storage size, hdr.dime.bitpix/8 will be
                  %  applied.
                  %
                  pos = sub2ind([d1 d2 d3 d4 d5 d6 d7], 1, 1, slice_idx(s), ...
			                    img_idx(t), dim5_idx(i5),dim6_idx(i6),dim7_idx(i7)) -1;
                  pos = pos * hdr.dime.bitpix/8;

                  % ROMAN: begin
                  if(roman)
                      % do nothing
                  else
                     img_siz = prod(hdr.dime.dim(2:3));

                     %  For complex float32 or complex float64, voxel values
                     %  include [real, imag]
                     %
                     if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
                        img_siz = img_siz * 2;
                     end

                     %MPH: For RGB24, voxel values include 3 separate color planes
                     %
                     if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
                        img_siz = img_siz * 3;
                     end
                  end; % if (roman)
                  % ROMAN: end
         
                  if filetype == 2
                     fseek(fid, pos + hdr.dime.vox_offset, 'bof');
                  else
                     fseek(fid, pos, 'bof');
                  end

                  %  For each frame, fread will read precision of value
                  %  in img_siz times
                  %
                  % ROMAN: begin
                  if(roman)
                     img(:,currentIndex) = fread(fid, img_siz, sprintf('*%s',precision));
                     currentIndex = currentIndex +1;
                  else
                     img = [img fread(fid, img_siz, sprintf('*%s',precision))];
                  end; %if(roman)
                  % ROMAN: end
                  
               end
               end
            end
         end
      end
   end
   
   %  For complex float32 or complex float64, voxel values
   %  include [real, imag]
   %
   if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
      img = reshape(img, [2, length(img)/2]);
      img = complex(img(1,:)', img(2,:)');
   end

   fclose(fid);

   %  Update the global min and max values 
   %
   hdr.dime.glmax = double(max(img(:)));
   hdr.dime.glmin = double(min(img(:)));

   %  old_RGB treat RGB slice by slice, now it is treated voxel by voxel
   %
   if old_RGB & hdr.dime.datatype == 128 & hdr.dime.bitpix == 24
      % remove squeeze
      img = (reshape(img, [hdr.dime.dim(2:3) 3 length(slice_idx) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
      img = permute(img, [1 2 4 3 5 6 7 8]);
   elseif hdr.dime.datatype == 128 & hdr.dime.bitpix == 24
      % remove squeeze
      img = (reshape(img, [3 hdr.dime.dim(2:3) length(slice_idx) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
      img = permute(img, [2 3 4 1 5 6 7 8]);
   elseif hdr.dime.datatype == 511 & hdr.dime.bitpix == 96
      img = double(img(:));
      img = single((img - min(img))/(max(img) - min(img)));
      % remove squeeze
      img = (reshape(img, [3 hdr.dime.dim(2:3) length(slice_idx) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
      img = permute(img, [2 3 4 1 5 6 7 8]);
   else
      % remove squeeze
      img = (reshape(img, [hdr.dime.dim(2:3) length(slice_idx) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
   end

   if ~isempty(slice_idx)
      hdr.dime.dim(4) = length(slice_idx);
   end

   if ~isempty(img_idx)
      hdr.dime.dim(5) = length(img_idx);
   end

   if ~isempty(dim5_idx)
      hdr.dime.dim(6) = length(dim5_idx);
   end

   if ~isempty(dim6_idx)
      hdr.dime.dim(7) = length(dim6_idx);
   end

   if ~isempty(dim7_idx)
      hdr.dime.dim(8) = length(dim7_idx);
   end

   return						% read_image

%------------------------------------------------
function save_untouch_nii(nii, filename)
   
   if ~exist('nii','var') | isempty(nii) | ~isfield(nii,'hdr') | ...
	~isfield(nii,'img') | ~exist('filename','var') | isempty(filename)

      error('Usage: save_untouch_nii(nii, filename)');
   end

   if ~isfield(nii,'untouch') | nii.untouch == 0
      error('Usage: please use ''save_nii.m'' for the modified structure.');
   end

   if isfield(nii.hdr.hist,'magic') & strcmp(nii.hdr.hist.magic(1:3),'ni1')
      filetype = 1;
   elseif isfield(nii.hdr.hist,'magic') & strcmp(nii.hdr.hist.magic(1:3),'n+1')
      filetype = 2;
   else
      filetype = 0;
   end

   [p,f] = fileparts(filename);
   fileprefix = fullfile(p, f);

   write_nii(nii, filetype, fileprefix);

%   %  So earlier versions of SPM can also open it with correct originator
 %  %
  % if filetype == 0
   %   M=[[diag(nii.hdr.dime.pixdim(2:4)) -[nii.hdr.hist.originator(1:3).*nii.hdr.dime.pixdim(2:4)]'];[0 0 0 1]];
    %  save(fileprefix, 'M');
%   elseif filetype == 1
 %     M=[];
  %    save(fileprefix, 'M');
   %end
   
   return					% save_untouch_nii


%-----------------------------------------------------------------------------------
function write_nii(nii, filetype, fileprefix)

   hdr = nii.hdr;

   if isfield(nii,'ext') & ~isempty(nii.ext)
      ext = nii.ext;
      [ext, esize_total] = verify_nii_ext(ext);
   else
      ext = [];
   end

   switch double(hdr.dime.datatype),
   case   1,
      hdr.dime.bitpix = int16(1 ); precision = 'ubit1';
   case   2,
      hdr.dime.bitpix = int16(8 ); precision = 'uint8';
   case   4,
      hdr.dime.bitpix = int16(16); precision = 'int16';
   case   8,
      hdr.dime.bitpix = int16(32); precision = 'int32';
   case  16,
      hdr.dime.bitpix = int16(32); precision = 'float32';
   case  32,
      hdr.dime.bitpix = int16(64); precision = 'float32';
   case  64,
      hdr.dime.bitpix = int16(64); precision = 'float64';
   case 128,
      hdr.dime.bitpix = int16(24); precision = 'uint8';
   case 256 
      hdr.dime.bitpix = int16(8 ); precision = 'int8';
   case 512 
      hdr.dime.bitpix = int16(16); precision = 'uint16';
   case 768 
      hdr.dime.bitpix = int16(32); precision = 'uint32';
   case 1024
      hdr.dime.bitpix = int16(64); precision = 'int64';
   case 1280
      hdr.dime.bitpix = int16(64); precision = 'uint64';
   case 1792,
      hdr.dime.bitpix = int16(128); precision = 'float64';
   otherwise
      error('This datatype is not supported');
   end
   
%   hdr.dime.glmax = round(double(max(nii.img(:))));
 %  hdr.dime.glmin = round(double(min(nii.img(:))));
   
   if filetype == 2
      fid = fopen(sprintf('%s.nii',fileprefix),'w');
      
      if fid < 0,
         msg = sprintf('Cannot open file %s.nii.',fileprefix);
         error(msg);
      end
      
      hdr.dime.vox_offset = 352;

      if ~isempty(ext)
         hdr.dime.vox_offset = hdr.dime.vox_offset + esize_total;
      end

      hdr.hist.magic = 'n+1';
      save_untouch_nii_hdr(hdr, fid);

      if ~isempty(ext)
         save_nii_ext(ext, fid);
      end
   elseif filetype == 1
      fid = fopen(sprintf('%s.hdr',fileprefix),'w');
      
      if fid < 0,
         msg = sprintf('Cannot open file %s.hdr.',fileprefix);
         error(msg);
      end
      
      hdr.dime.vox_offset = 0;
      hdr.hist.magic = 'ni1';
      save_untouch_nii_hdr(hdr, fid);

      if ~isempty(ext)
         save_nii_ext(ext, fid);
      end
      
      fclose(fid);
      fid = fopen(sprintf('%s.img',fileprefix),'w');
   else
      fid = fopen(sprintf('%s.hdr',fileprefix),'w');
      
      if fid < 0,
         msg = sprintf('Cannot open file %s.hdr.',fileprefix);
         error(msg);
      end
      
      save_untouch0_nii_hdr(hdr, fid);
      
      fclose(fid);
      fid = fopen(sprintf('%s.img',fileprefix),'w');
   end

   ScanDim = double(hdr.dime.dim(5));		% t
   SliceDim = double(hdr.dime.dim(4));		% z
   RowDim   = double(hdr.dime.dim(3));		% y
   PixelDim = double(hdr.dime.dim(2));		% x
   SliceSz  = double(hdr.dime.pixdim(4));
   RowSz    = double(hdr.dime.pixdim(3));
   PixelSz  = double(hdr.dime.pixdim(2));
   
   x = 1:PixelDim;
   
   if filetype == 2 & isempty(ext)
      skip_bytes = double(hdr.dime.vox_offset) - 348;
   else
      skip_bytes = 0;
   end

   if double(hdr.dime.datatype) == 128

      %  RGB planes are expected to be in the 4th dimension of nii.img
      %
      if(size(nii.img,4)~=3)
         error(['The NII structure does not appear to have 3 RGB color planes in the 4th dimension']);
      end

      nii.img = permute(nii.img, [4 1 2 3 5 6 7 8]);
   end

   %  For complex float32 or complex float64, voxel values
   %  include [real, imag]
   %
   if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
      real_img = real(nii.img(:))';
      nii.img = imag(nii.img(:))';
      nii.img = [real_img; nii.img];
   end

   if skip_bytes
      fwrite(fid, zeros(1,skip_bytes), 'uint8');
   end

   fwrite(fid, nii.img, precision);
%   fwrite(fid, nii.img, precision, skip_bytes);        % error using skip
   fclose(fid);

   return;					% write_nii




%------------------------------------------------
function save_untouch_nii_hdr(hdr, fid)

   if ~isequal(hdr.hk.sizeof_hdr,348),
      error('hdr.hk.sizeof_hdr must be 348.');
   end

   write_header(hdr, fid);

   return;					% save_nii_hdr


%---------------------------------------------------------------------
function write_header(hdr, fid)

        %  Original header structures
	%  struct dsr				/* dsr = hdr */
	%       { 
	%       struct header_key hk;            /*   0 +  40       */
	%       struct image_dimension dime;     /*  40 + 108       */
	%       struct data_history hist;        /* 148 + 200       */
	%       };                               /* total= 348 bytes*/
   
   header_key2(fid, hdr.hk);
   image_dimension2(fid, hdr.dime);
   data_history2(fid, hdr.hist);
   
   %  check the file size is 348 bytes
   %
   fbytes = ftell(fid);
   
   if ~isequal(fbytes,348),
      msg = sprintf('Header size is not 348 bytes.');
      warning(msg);
   end
    
   return;					% write_header


%---------------------------------------------------------------------
function header_key2(fid, hk)
   
   fseek(fid,0,'bof');

	%  Original header structures    
	%  struct header_key                      /* header key      */ 
	%       {                                /* off + size      */
	%       int sizeof_hdr                   /*  0 +  4         */
	%       char data_type[10];              /*  4 + 10         */
	%       char db_name[18];                /* 14 + 18         */
	%       int extents;                     /* 32 +  4         */
	%       short int session_error;         /* 36 +  2         */
	%       char regular;                    /* 38 +  1         */
	%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
	%       };                               /* total=40 bytes  */
        
   fwrite(fid, hk.sizeof_hdr(1),    'int32');	% must be 348.
    
   % data_type = sprintf('%-10s',hk.data_type);	% ensure it is 10 chars from left
   % fwrite(fid, data_type(1:10), 'uchar');
   pad = zeros(1, 10-length(hk.data_type));
   hk.data_type = [hk.data_type  char(pad)];
   fwrite(fid, hk.data_type(1:10), 'uchar');
    
   % db_name   = sprintf('%-18s', hk.db_name);	% ensure it is 18 chars from left
   % fwrite(fid, db_name(1:18), 'uchar');
   pad = zeros(1, 18-length(hk.db_name));
   hk.db_name = [hk.db_name  char(pad)];
   fwrite(fid, hk.db_name(1:18), 'uchar');
    
   fwrite(fid, hk.extents(1),       'int32');
   fwrite(fid, hk.session_error(1), 'int16');
   fwrite(fid, hk.regular(1),       'uchar');	% might be uint8
    
   % fwrite(fid, hk.hkey_un0(1),    'uchar');
   % fwrite(fid, hk.hkey_un0(1),    'uint8');
   fwrite(fid, hk.dim_info(1),      'uchar');
    
   return;					% header_key


%---------------------------------------------------------------------
function image_dimension2(fid, dime)

	%  Original header structures        
	%  struct image_dimension
	%       {                                /* off + size      */
	%       short int dim[8];                /* 0 + 16          */
	%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
	%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
	%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
	%       short int intent_code;   % short int unused1;   /* 28 + 2 */
	%       short int datatype;              /* 30 + 2          */
	%       short int bitpix;                /* 32 + 2          */
	%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
	%       float pixdim[8];                 /* 36 + 32         */
	%			/*
	%				pixdim[] specifies the voxel dimensions:
	%				pixdim[1] - voxel width
	%				pixdim[2] - voxel height
	%				pixdim[3] - interslice distance
	%				pixdim[4] - volume timing, in msec
	%					..etc
	%			*/
	%       float vox_offset;                /* 68 + 4          */
	%       float scl_slope;   % float roi_scale;     /* 72 + 4 */
	%       float scl_inter;   % float funused1;      /* 76 + 4 */
	%       short slice_end;   % float funused2;      /* 80 + 2 */
	%       char slice_code;   % float funused2;      /* 82 + 1 */
	%       char xyzt_units;   % float funused2;      /* 83 + 1 */
	%       float cal_max;                   /* 84 + 4          */
	%       float cal_min;                   /* 88 + 4          */
	%       float slice_duration;   % int compressed; /* 92 + 4 */
	%       float toffset;   % int verified;          /* 96 + 4 */
	%       int glmax;                       /* 100 + 4         */
	%       int glmin;                       /* 104 + 4         */
	%       };                               /* total=108 bytes */
	
   fwrite(fid, dime.dim(1:8),        'int16');
   fwrite(fid, dime.intent_p1(1),  'float32');
   fwrite(fid, dime.intent_p2(1),  'float32');
   fwrite(fid, dime.intent_p3(1),  'float32');
   fwrite(fid, dime.intent_code(1),  'int16');
   fwrite(fid, dime.datatype(1),     'int16');
   fwrite(fid, dime.bitpix(1),       'int16');
   fwrite(fid, dime.slice_start(1),  'int16');
   fwrite(fid, dime.pixdim(1:8),   'float32');
   fwrite(fid, dime.vox_offset(1), 'float32');
   fwrite(fid, dime.scl_slope(1),  'float32');
   fwrite(fid, dime.scl_inter(1),  'float32');
   fwrite(fid, dime.slice_end(1),    'int16');
   fwrite(fid, dime.slice_code(1),   'uchar');
   fwrite(fid, dime.xyzt_units(1),   'uchar');
   fwrite(fid, dime.cal_max(1),    'float32');
   fwrite(fid, dime.cal_min(1),    'float32');
   fwrite(fid, dime.slice_duration(1), 'float32');
   fwrite(fid, dime.toffset(1),    'float32');
   fwrite(fid, dime.glmax(1),        'int32');
   fwrite(fid, dime.glmin(1),        'int32');
   
   return;					% image_dimension


%---------------------------------------------------------------------
function data_history2(fid, hist)
    
	% Original header structures
	%struct data_history       
	%       {                                /* off + size      */
	%       char descrip[80];                /* 0 + 80          */
	%       char aux_file[24];               /* 80 + 24         */
	%       short int qform_code;            /* 104 + 2         */
	%       short int sform_code;            /* 106 + 2         */
	%       float quatern_b;                 /* 108 + 4         */
	%       float quatern_c;                 /* 112 + 4         */
	%       float quatern_d;                 /* 116 + 4         */
	%       float qoffset_x;                 /* 120 + 4         */
	%       float qoffset_y;                 /* 124 + 4         */
	%       float qoffset_z;                 /* 128 + 4         */
	%       float srow_x[4];                 /* 132 + 16        */
	%       float srow_y[4];                 /* 148 + 16        */
	%       float srow_z[4];                 /* 164 + 16        */
	%       char intent_name[16];            /* 180 + 16        */
	%       char magic[4];   % int smin;     /* 196 + 4         */
	%       };                               /* total=200 bytes */
	
   % descrip     = sprintf('%-80s', hist.descrip);     % 80 chars from left
   % fwrite(fid, descrip(1:80),    'uchar');
   pad = zeros(1, 80-length(hist.descrip));
   hist.descrip = [hist.descrip  char(pad)];
   fwrite(fid, hist.descrip(1:80), 'uchar');
    
   % aux_file    = sprintf('%-24s', hist.aux_file);    % 24 chars from left
   % fwrite(fid, aux_file(1:24),   'uchar');
   pad = zeros(1, 24-length(hist.aux_file));
   hist.aux_file = [hist.aux_file  char(pad)];
   fwrite(fid, hist.aux_file(1:24), 'uchar');
    
   fwrite(fid, hist.qform_code,    'int16');
   fwrite(fid, hist.sform_code,    'int16');
   fwrite(fid, hist.quatern_b,   'float32');
   fwrite(fid, hist.quatern_c,   'float32');
   fwrite(fid, hist.quatern_d,   'float32');
   fwrite(fid, hist.qoffset_x,   'float32');
   fwrite(fid, hist.qoffset_y,   'float32');
   fwrite(fid, hist.qoffset_z,   'float32');
   fwrite(fid, hist.srow_x(1:4), 'float32');
   fwrite(fid, hist.srow_y(1:4), 'float32');
   fwrite(fid, hist.srow_z(1:4), 'float32');

   % intent_name = sprintf('%-16s', hist.intent_name);	% 16 chars from left
   % fwrite(fid, intent_name(1:16),    'uchar');
   pad = zeros(1, 16-length(hist.intent_name));
   hist.intent_name = [hist.intent_name  char(pad)];
   fwrite(fid, hist.intent_name(1:16), 'uchar');
    
   % magic	= sprintf('%-4s', hist.magic);		% 4 chars from left
   % fwrite(fid, magic(1:4),           'uchar');
   pad = zeros(1, 4-length(hist.magic));
   hist.magic = [hist.magic  char(pad)];
   fwrite(fid, hist.magic(1:4),        'uchar');
    
   return;					% data_history

%------------------------------------------------
function [ext, esize_total] = verify_nii_ext(ext)

   if ~isfield(ext, 'section')
      error('Incorrect NIFTI header extension structure.');
   elseif ~isfield(ext, 'num_ext')
      ext.num_ext = length(ext.section);
   elseif ~isfield(ext, 'extension')
      ext.extension = [1 0 0 0];
   end

   esize_total = 0;

   for i=1:ext.num_ext
      if ~isfield(ext.section(i), 'ecode') | ~isfield(ext.section(i), 'edata')
         error('Incorrect NIFTI header extension structure.');
      end

      ext.section(i).esize = ceil((length(ext.section(i).edata)+8)/16)*16;
      ext.section(i).edata = ...
	[ext.section(i).edata ...
	 zeros(1,ext.section(i).esize-length(ext.section(i).edata)-8)];
      esize_total = esize_total + ext.section(i).esize;
   end

   return                                       % verify_nii_ext

%------------------------------------------------
function save_nii_ext(ext, fid)

   if ~exist('ext','var') | ~exist('fid','var')
      error('Usage: save_nii_ext(ext, fid)');
   end

   if ~isfield(ext,'extension') | ~isfield(ext,'section') | ~isfield(ext,'num_ext')
      error('Wrong header extension');
   end

   write_ext(ext, fid);

   return;                                      % save_nii_ext


%---------------------------------------------------------------------
function write_ext(ext, fid)

   fwrite(fid, ext.extension, 'uchar');

   for i=1:ext.num_ext
      fwrite(fid, ext.section(i).esize, 'int32');
      fwrite(fid, ext.section(i).ecode, 'int32');
      fwrite(fid, ext.section(i).edata, 'uchar');
   end

   return;                                      % write_ext

%------------------------------------------------








