% Usage:cbarref = overlay_jet(baseimage,map,mask,map_range,baseimagemax(optional));
%
% NOTE: map_range must be symmetrical about 0 for red = +ve, blue = -ve

function cbarref = overlay_jet(baseimage,map,mask,map_range,baseimagemax)

cm = colormap(jet(128));

%% Removing NaN values from map (and masking out any NaN values in plot)
map(isnan(map))=0;
mask(isnan(map))=0;

%% Scaling base image between 0 and 1, for true colour format
if nargin == 5
    baseimage = (baseimage-min(baseimage(:)))./(baseimagemax-min(baseimage(:)));
else
    baseimage = (baseimage-min(baseimage(:)))./(max(baseimage(:))-min(baseimage(:)));
end
baseimage(baseimage<0)=0;baseimage(baseimage>1)=1;

%% Rotating images

s = size(baseimage);
if size(s,2) < 3
	s(3) = 1;
end

rbaseimage = zeros(s(2),s(1),s(3));
for n = 1:s(3)
	rbaseimage(:,:,n) = baseimage(:,:,n)';% swapping x and y coordinates
end

rbaseimage = rbaseimage((1:end)*-1+s(2)+1,:,:);% Inverts x coordinate
rbaseimage = reshape(rbaseimage,s(2),s(1)*s(3));


rmap = zeros(s(2),s(1),s(3));
for n = 1:s(3)
	rmap(:,:,n) = map(:,:,n)';
end

rmap = rmap((1:end)*-1+s(2)+1,:,:);
rmap = reshape(rmap,s(2),s(1)*s(3));

rmask = zeros(s(2),s(1),s(3));
for n = 1:s(3)
	rmask(:,:,n) = mask(:,:,n)';
end

rmask = rmask((1:end)*-1+s(2)+1,:,:);
rmask = reshape(rmask,s(2),s(1)*s(3));

%% 3 column, truecolour (rgb) format

if nargin > 3 && numel(map_range)~=0
	rmap(rmap<map_range(1))=map_range(1);rmap(rmap>map_range(2))=map_range(2);
else
	map_range = [min(rmap(rmask~=0)) max(rmap(rmask~=0))];
end

% Assigning each map value to a 3 element rgb true colour.
img = cm(round(127*(rmap-map_range(1))./(map_range(2)-map_range(1))+1),:);
img = reshape(img,s(2),s(1)*s(3),3);

img(repmat(rmask,[1 1 3])==0) = repmat(rbaseimage(rmask==0),[1 1 3]);% Pixels outside mask show baseimage in greyscale

imagesc(img)

cbarref = [linspace(map_range(1),map_range(2),128)' cm];
