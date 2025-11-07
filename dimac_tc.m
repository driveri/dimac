function out = dimac_tc(dimac_fname,outprefix_mask)
% Function to generate DIMAC timecourse data.
% IDD 16/10/2025
%
% Accepts DIMAC 4D data; displays DIMAC timeseries so the user can choose the appropriate ROI for the output timecourse.
%
% Usage: out = dimac_tc(dimac_fname);
%
%        out          - output structure including mean timecourse over ROI (out.tc), TR of DIMAC data (out.tr)
%
%        dimac_fname  - input DIMAC dataset filename (include full extension of nifti file)
%
%      outprefix_mask - [OPTIONAL INPUT] filename prefix string to save resulting mask as nifti (added 25/09/2024 IDD).


%% Load in DIMAC 4D data:
dimac = load_untouch_nii(dimac_fname);

%% Initialise cluster map and output:
cm = ones(size(dimac.img,1),size(dimac.img,2),size(dimac.img,3));
out = [];

tcs = [];
grad_test = [];

pc1 = 90; % Initial percentile used to set signal intensity threshold for the ROIs
c_ind = 0;

% IDD 25/09/2024 - added optional third input for defining filename for saving mask as nifti
if nargin > 1
    % Third input indicates saving mask as nifti
    omaskflag=true;
else
    omaskflag=false;
end

%% Initialising GUI:
f1 = figure('Position',[100 50 1000 450]);
a1 = axes('units','pixels','Position',[25 25 400 400]);
a2 = axes('units','pixels','Position',[475 25 500 200],'visible','off');
a3 = axes('units','pixels','Position',[675 275 150 150],'visible','off');
h_percentile_input = uicontrol('Style','edit','Position',[850 350 50 20],...
    'String',pc1,'Callback', @plotROItcs);
%h_popupchooseROI = uicontrol('Style', 'popup', 'String', 0 ,...
%    'Position', [475 350 100 50], 'visible', 'on',...
%    'Callback', @plotROItcs);
h_chooseROI = uicontrol('Style', 'togglebutton', 'String', 'Choose ROI' ,...
    'Position', [475 350 100 50], 'visible', 'on',...
    'Callback', @plotROItcs);
h_roisel_txt = uicontrol('Style','text',...
    'Position',[450 400 150 20],...
    'String','Choose artery','visible','off');
h_out_btn = uicontrol('Style', 'pushbutton', 'String', 'Continue and Close',...
    'Position', [475 300 150 50],'visible','off',...
    'Callback', @buttonfunction);
h_plot_xrange = uicontrol('Style','edit','Position',[550 230 50 20],...
    'visible','off','String',size(dimac.img,4),'Callback', @plotxlim);
h_plot_xmin = uicontrol('Style','edit','Position',[475 230 50 20],...
    'visible','off','String',0,'Callback', @plotxlim);

h_bg = uibuttongroup(f1,'Visible','off',...
    'Position',[0.58 0.83 0.08 0.12],...
    'SelectionChangedFcn',@plotROItcs);
% [525 350 100 100]
% Create three radio buttons in the button group.
h_radio_peak25 = uicontrol(h_bg,'Style','radiobutton',...
    'String','peak',...
    'Position',[10 5 80 20],...
    'HandleVisibility','off');
h_radio_cluster = uicontrol(h_bg,'Style',...
    'radiobutton',...
    'String','cluster',...
    'Position',[10 30 80 20],...
    'HandleVisibility','off');

h_bg.Visible = 'on';

h_gradtest_tickbox = uicontrol('Style','checkbox','Position',[850 400 20 20],...
    'visible','on','Callback', @plotROItcs);
h_gradtest_plotbutton = uicontrol('Style','pushbutton',...
    'string','gradient test','Position',[880 400 100 20],...
    'visible','on','Callback',@plotgradfit);

if nargin > 1
    if numel(outprefix_mask)>41
        annotation('textbox',[0.03 0.88 0.2 0.05],'String',regexprep(outprefix_mask(end-41:end),'_',' '),'Color',[1 1 1],'LineStyle','none','FontSize',12)
    else
        annotation('textbox',[0.03 0.88 0.2 0.05],'String',regexprep(outprefix_mask,'_',' '),'Color',[1 1 1],'LineStyle','none','FontSize',12)
    end
end

% Calculate cluster map based on pc1 percentile threshold and watershed
% segmentation:
cm = cluster_on_percentile(mean(dimac.img,4),pc1);
mask1 = cm>0;
% h_popupchooseROI.String = unique(cm(cm~=0))';


waitfor(f1) % wait for figure to close before returning function output

    function plotROItcs(source,callbackdata)

        
        if sum(mask1(:))==0
            axes(a2);cla;a2.Visible = 'off';
            set(h_out_btn,'visible','off')
            axes(a3);overlay_mask(mean(dimac.img,4));axis off
        else
            % Check which ui is calling the function
            if h_chooseROI.Value == 1 % strcmp(class(source),'matlab.ui.control.UIControl') && strcmp(source.Style,'togglebutton')
                %disp('roi selection')
                poi = drawpoint(a1); % Select the ROI of interest - will use these coordinates to search for the nearest ROI
                exitflag = false;
                while ~exitflag
                    search_rad = 0;
                    c_ind = nonzeros(cm(round(poi.Position(1))-search_rad:round(poi.Position(1))+search_rad,round(100-poi.Position(2))-search_rad:round(100-poi.Position(2))+search_rad));
                    while numel(c_ind)==0
                        search_rad = search_rad + 1; % Expand search window until finds an ROI (cm)
                        c_ind = nonzeros(cm(round(poi.Position(1))-search_rad:round(poi.Position(1))+search_rad,round(100-poi.Position(2))-search_rad:round(100-poi.Position(2))+search_rad));
                    end
                    if numel(unique(c_ind))==1 % Check in case two ROIs are selected - if so, will wait until a new position selected, where only a single ROI is resolved
                        exitflag = true;
                    end
                    pause(0.1)
                    
                end
                poi.Visible = 'off';
                %c_ind = str2double(source.String(source.Value,:));
                h_chooseROI.Value = 0;
            else
                % Update cluster list if not called by the list
                pc1 = str2double(h_percentile_input.String);
                cm = cluster_on_percentile(mean(dimac.img,4),pc1);
                mask1 = cm>0;
%                h_popupchooseROI.String = unique(cm(cm~=0))';
%                if c_ind > max(cm(:)) % If number of clusters reduced, need to reset the current value to be within the new range, so the list is rendered
%                    h_popupchooseROI.Value = 1;
%                end
            end
%            if c_ind ~=0 % only apply the cluster if index ~= 0 (which means plot whole mask)
                mask1 = cm==c_ind;
                if h_radio_peak25.Value == 1
                    mask1 = imdilate(mask1,ones(5,5));
                    % If using the peak voxel, choose a 3x3 ROI, centred on
                    % the selected peak voxel.
                end
%            end

            
            %% Plotting ROI overlay and timecourse for selected ROI:
            %     mask1 = km.img==km_ind;
            tcs = reshape(dimac.img(repmat(mask1,[1 1 1 size(dimac.img,4)])),[],size(dimac.img,4))';
            tr = dimac.hdr.dime.pixdim(5); % TR in seconds, read from the nifti header
            if h_gradtest_tickbox.Value==1 % if choose to perform gradient test (see dimac_gradientsearch.m)
                grad_test = dimac_gradientsearch(tcs,tr);
                mask1(nonzeros(find(mask1).*~grad_test.mask)) = false; % Removing voxels from mask1 which don't pass the gradient test (final slope within 0.5-1.5x the slope of the total slope for > half of beats)
                
            end
            if sum(mask1(:))>0 % only calculate and plot timeseries if the mask is not empty
                tc = mean(reshape(dimac.img(repmat(mask1,[1 1 1 size(dimac.img,4)])),[],size(dimac.img,4)),1)';
                ppr = pulsepower(tc,tr); % Pulse power ratio (ratio of power in 40-120bpm range to total power)
                axes(a3);overlay_jet(mean(dimac.img,4),cm,mask1,[0 max(cm(:))/3]);axis off
                axes(a2);plot(reshape(dimac.img(repmat(mask1,[1 1 1 size(dimac.img,4)])),[],size(dimac.img,4))');
                hold on;plot(tc,'k-','LineWidth',2);hold off;title(['Pulse power ratio = ',num2str(ppr)])
                box off
                plotxlim
                a3.Visible = 'on';a3.XTickLabel = {};a3.YTickLabel = {};
                
                h_roisel_txt.Visible = 'on';
                h_popupchooseROI.Visible = 'on';
                if c_ind ~=0
                    h_plot_xmin.Visible = 'on';
                    h_plot_xrange.Visible = 'on';
                    set(h_out_btn,'visible','on')
                end
                out.tc = tc;
                out.tr = tr;
                out.ppr = ppr;
                out.mask = mask1;
                
            end
            
        end
    end

    function cm = cluster_on_percentile(im3D,pc1)
        if h_radio_cluster.Value == 1
            cm = double(watershed(-im3D.*(im3D>prctile(reshape(im3D,numel(im3D),[]),pc1)))).*(im3D>prctile(reshape(im3D,numel(im3D),[]),pc1));
        elseif h_radio_peak25.Value == 1
            
            peaks1 = imregionalmax(im3D.*(im3D>prctile(reshape(im3D,numel(im3D),[]),pc1)),8);
            cm = double(watershed(-peaks1)).*peaks1;clear peaks1
        else
            error('Unexpected option from the button group h_bg')
        end
        % Initial plot
        axes(a1);overlay_jet(mean(dimac.img,4),cm,cm>0,[0 max(cm(:))]);
        axis off 
    end

    function buttonfunction(source,callbackdata)
%        [out.peak,out.base,out.fit2,out.coeffs,out.R2,out.footind]=dimac_peak_extract(out.tc,numel(out.tc),out.tr);
%        out.K = 5; % Currently hard-coded into dimac_peak_extract.m; K is the # of Fourier pairs fit to each beat.
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

    function plotgradfit(source,callbackdata)
        if h_gradtest_tickbox.Value==1
            figure('WindowState','maximized');
            for n = 1:25
                subplot(5,5,n)
                plot(tcs(:,n))
                xlim([0 500])
                hold on
                for beatnum = 1:7
                    plot([grad_test.maxind(beatnum) grad_test.maxind(beatnum+1)],[grad_test.maxind(beatnum) grad_test.maxind(beatnum+1)]*grad_test.gradfit(beatnum,n,2)+grad_test.gradfit(beatnum,n,1),'r')
                    plot([grad_test.maxind(beatnum) grad_test.maxind(beatnum+1)],[grad_test.maxind(beatnum) grad_test.maxind(beatnum+1)]*grad_test.gradallfit(beatnum,n,2)+grad_test.gradallfit(beatnum,n,1),'g')
                end
                hold off
                if grad_test.mask(n)
                    title([num2str(mean(grad_test.gradratio(:,n))),' +/- ',num2str(std(grad_test.gradratio(:,n))),'; ',num2str(sum((grad_test.gradratio(:,n)>0.4).*(grad_test.gradratio(:,n)<2).*(grad_test.gradfit(:,n,2)<0),1)),'/',num2str(numel(grad_test.maxind)-1),' beats'],'Color','blue')
                else
                    title([num2str(mean(grad_test.gradratio(:,n))),' +/- ',num2str(std(grad_test.gradratio(:,n))),'; ',num2str(sum((grad_test.gradratio(:,n)>0.4).*(grad_test.gradratio(:,n)<2).*(grad_test.gradfit(:,n,2)<0),1)),'/',num2str(numel(grad_test.maxind)-1),' beats'])
                end
            end
        end
    end
end % Nested 2nd-level functions inside top function, so can pass variables down to them (N.B. cannot pass variables directly between 2nd-level functions)
