function out = dimac_tc(dimac_fname,kmeans_fname,outprefix_mask)
% Function to generate DIMAC timecourse data.
% IDD 13/08/2024
%
% Accepts DIMAC 4D data and indexed cluster masks (e.g from 3dkmeans); displays DIMAC timeseries so the user can choose the appropriate ROI for the output timecourse.
%
% Usage: out = dimac_tc(dimac_fname,kmeans_fname);
%
%        out          - output structure including mean timecourse over ROI (out.tc), TR of DIMAC data (out.tr)
%
%        dimac_fname  - input DIMAC dataset filename (include full extension of nifti file)
%        kmeans_fname - input indexed cluster mask filename (include full extension of nifti file)
%                          - script will discriminate clusters based on integer values
%
%      outprefix_mask - [OPTIONAL INPUT] filename prefix string to save resulting mask as nifti (added 25/09/2024 IDD).


% addpath(genpath('/home/sapid1/DIMAC_scripts'))
% Pointing to load_untouch_nii.m, pulsepower.m

% addpath ~/matlab/overlay_scripts/
% Pointing to image overlay plotting scripts: overlay_mask.m
% IDD 10/09/2024 - made a local copy of overlay_mask.m, so no longer need this local path

%% Load in DIMAC 4D data and k-means mask:
dimac = load_untouch_nii(dimac_fname);
km = load_untouch_nii(kmeans_fname);

% map to keep track of current cluster (to determine whether to edit
% current cluster, or calculate new clusters (based on order of operations)
current_cluster = zeros(size(km.img));
cm = ones(size(km.img));

%% Initialise output:
out = [];

% IDD 25/09/2024 - added optional third input for defining filename for saving mask as nifti
if nargin > 2
    % Third input indicates saving mask as nifti
    omaskflag=true;
else
    omaskflag=false;
end

%% Initialising GUI:
f1 = figure('Position',[100 50 1000 450]);
a1 = axes('units','pixels','Position',[25 25 400 400]);
a2 = axes('units','pixels','Position',[475 25 500 200],'visible','off');
a3 = axes('units','pixels','Position',[675 275 275 150]);
h_popupchooseROI = uicontrol('Style', 'popup', 'String', 0 ,...
    'Position', [475 350 100 50], 'visible', 'off',...
    'Callback', @plotROItcs);
h_roisel_txt = uicontrol('Style','text',...
    'Position',[475 400 150 20],...
    'String','Choose artery (0 = all)','visible','off');
h_out_btn = uicontrol('Style', 'pushbutton', 'String', 'Continue and Close',...
    'Position', [475 300 150 50],'visible','off',...
    'Callback', @buttonfunction);
h_plot_xrange = uicontrol('Style','edit','Position',[550 230 50 20],...
    'visible','off','String',size(dimac.img,4),'Callback', @plotxlim);
h_plot_xmin = uicontrol('Style','edit','Position',[475 230 50 20],...
    'visible','off','String',0,'Callback', @plotxlim);

for init_count = 1:floor(max(km.img(:)))
    h_checkbox_kmeans(init_count) = uicontrol('Style','checkbox','String',num2str(init_count),'Position',[440 400-30*(init_count-1) 30 20],'Callback', @plotROItcs);
end
if nargin > 2
    if numel(outprefix_mask)>41
        annotation('textbox',[0.03 0.88 0.2 0.05],'String',regexprep(outprefix_mask(end-41:end),'_',' '),'Color',[1 1 1],'LineStyle','none','FontSize',12)
    else
        annotation('textbox',[0.03 0.88 0.2 0.05],'String',regexprep(outprefix_mask,'_',' '),'Color',[1 1 1],'LineStyle','none','FontSize',12)
    end
end

axes(a1)
overlay_mask(mean(dimac.img,4))
axis off
% Plot each ROI timeseries for overview:
tc_survey = zeros(size(dimac.img,4),floor(max(km.img(:))));
for roi_count = 1:floor(max(km.img(:)))
    mask1 = km.img==roi_count;
    tc_survey(:,roi_count) = mean(reshape(dimac.img(repmat(mask1,[1 1 1 size(dimac.img,4)])),[],size(dimac.img,4)),1)';
    clear mask1
end
axes(a3)
plot(tc_survey(1:500,:))
legend('Location','eastoutside')
axis off
waitfor(f1) % wait for figure to close before returning function output

    function plotROItcs(source,callbackdata)
        
        
        mask1 = zeros(size(km.img));
        for n = 1:floor(max(km.img(:)))
            mask1 = mask1+(double(km.img)==n)*h_checkbox_kmeans(n).Value;
        end
        mask1 = mask1>0; % convert to logical format, for use with indexing
        
        
        
        
        % km_ind = str2double(source.String(source.Value));
        if sum(mask1(:))==0
            axes(a2);cla;a2.Visible = 'off';
            set(h_out_btn,'visible','off')
            axes(a1);overlay_mask(mean(dimac.img,4))
        else
            % Check which ui is calling the function
            if strcmp(source.Style,'popupmenu')
                c_ind = str2double(source.String(source.Value,:));
                
            else
                c_ind = 0; % Update cluster list if not called by the list
                h_popupchooseROI.Value = 1;
                cm = connectedfun2D(mask1); % Simple connectivity analysis (searches with a 3x3 kernel)
                h_popupchooseROI.String = unique(cm(:))';
            end
            if c_ind ~=0 % only apply the cluster if index ~= 0 (which means plot whole mask)
                mask1 = cm==c_ind;
            end
            
            %% Plotting ROI overlay and timecourse for selected ROI:
            %     mask1 = km.img==km_ind;
            tc = mean(reshape(dimac.img(repmat(mask1,[1 1 1 size(dimac.img,4)])),[],size(dimac.img,4)),1)';
            tr = dimac.hdr.dime.pixdim(5); % TR in seconds, read from the nifti header
            ppr = pulsepower(tc,tr); % Pulse power ratio (ratio of power in 40-120bpm range to total power)
            axes(a1);overlay_jet(mean(dimac.img,4),cm,mask1,[0 max(cm(:))]);
            axes(a2);plot(reshape(dimac.img(repmat(mask1,[1 1 1 size(dimac.img,4)])),[],size(dimac.img,4))');
            hold on;plot(tc,'k-','LineWidth',2);hold off;title(['Pulse power ratio = ',num2str(ppr)])
            box off
            plotxlim
            h_roisel_txt.Visible = 'on';
            h_popupchooseROI.Visible = 'on';
            h_plot_xmin.Visible = 'on';
            h_plot_xrange.Visible = 'on';
            set(h_out_btn,'visible','on')
            out.tc = tc;
            out.tr = tr;
            out.ppr = ppr;
            out.mask = mask1;
        end
    end

    function buttonfunction(source,callbackdata)
        [out.peak,out.base,out.fit2,out.coeffs,out.R2,out.footind]=dimac_peak_extract(out.tc,numel(out.tc),out.tr);
        out.K = 5; % Currently hard-coded into dimac_peak_extract.m; K is the # of Fourier pairs fit to each beat.
        close(gcf)
        
        % IDD 25/09/2024 - added optional third input for defining filename for saving mask as nifti
        if omaskflag
            msk1 = dimac; % Copying header information from DIMAC image
            msk1.img = out.mask; % Replace image data with dilated (4x4) mask
            msk1.hdr.dime.dim([1 5]) = [3 1]; % reshape header to 3D
            save_untouch_nii(msk1,[outprefix_mask,'_roi.nii.gz'])
        end
        
    end

    function plotxlim(source,callbackdata)
        xmin = str2double(h_plot_xmin.String);
        xrange = str2double(h_plot_xrange.String);
        a2.XLim = [xmin xmin+xrange];
    end
end % Nested 2nd-level functions inside top function, so can pass variables down to them (N.B. cannot pass variables directly between 2nd-level functions)
