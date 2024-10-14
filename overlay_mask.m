% Function to display up to 5 binary masks on top of a greyscale image
%
% overlay_mask(img,mask1,mask2,mask3,mask4,mask5);
%
% mask1 = red; mask2 = blue; mask3 = green; mask4 = yellow; mask5 = cyan
%
%
% To Do:
% Replace looping with permute.m
% imagesc(reshape(permute(reshape(permute(reshape(img(:,64:-1:1,15:39),[64 64 5 5]),[1 2 4 3]),64,64*5,5),[2 1 3]),64*5,64*5))

function overlay_mask(img,mask1,mask2,mask3,mask4,mask5)

s = size(img);
if size(s,2) < 3
	s(3) = 1;
end
%% rotating image
for n = 1:s(3)
	rimg(:,:,n) = img(:,:,n)';% swapping x and y coordinates
end

rimg = rimg((1:end)*-1+s(2)+1,:,:);% Inverts x coordinate
rimg = reshape(rimg,s(2),s(1)*s(3));

im1 = imagesc(rimg);
colormap(gray)
hold on

if nargin > 1
	for n = 1:s(3)
		rmask1(:,:,n) = mask1(:,:,n)';
	end

	rmask1 = rmask1((1:end)*-1+s(2)+1,:,:);
	rmask1 = reshape(rmask1,s(2),s(1)*s(3));

	im2 = imagesc(rmask1);
	cmap1 = rmask1;cmap1(:,:,2:3) = zeros(s(2),s(1)*s(3),2);
	set(im2,'AlphaData',rmask1*0.5,'CData',cmap1)
end

if nargin > 2
	for n = 1:s(3)
		rmask2(:,:,n) = mask2(:,:,n)';
	end

	rmask2 = rmask2((1:end)*-1+s(2)+1,:,:);
	rmask2 = reshape(rmask2,s(2),s(1)*s(3));

	im3 = imagesc(rmask2);
	cmap2 = zeros(s(2),s(1)*s(3),3);cmap2(:,:,3) = rmask2;
	set(im3,'AlphaData',rmask2*0.5,'CData',cmap2)
end

if nargin > 3
	for n = 1:s(3)
		rmask3(:,:,n) = mask3(:,:,n)';
	end

	rmask3 = rmask3((1:end)*-1+s(2)+1,:,:);
	rmask3 = reshape(rmask3,s(2),s(1)*s(3));

	im4 = imagesc(rmask3);
	cmap3 = zeros(s(2),s(1)*s(3),3);cmap3(:,:,2) = rmask3;
	set(im4,'AlphaData',rmask3*0.5,'CData',cmap3)
end

if nargin > 4
	for n = 1:s(3)
		rmask4(:,:,n) = mask4(:,:,n)';
	end

	rmask4 = rmask4((1:end)*-1+s(2)+1,:,:);
	rmask4 = reshape(rmask4,s(2),s(1)*s(3));

	im5 = imagesc(rmask4);
	cmap4 = zeros(s(2),s(1)*s(3),3);cmap4(:,:,1) = rmask4;cmap4(:,:,2) = rmask4;
	set(im5,'AlphaData',rmask4*0.5,'CData',cmap4)
end

if nargin > 5
	for n = 1:s(3)
		rmask5(:,:,n) = mask5(:,:,n)';
	end

	rmask5 = rmask5((1:end)*-1+s(2)+1,:,:);
	rmask5 = reshape(rmask5,s(2),s(1)*s(3));

	im6 = imagesc(rmask5);
	cmap5 = zeros(s(2),s(1)*s(3),3);cmap5(:,:,2) = rmask5;cmap5(:,:,3) = rmask5;
	set(im6,'AlphaData',rmask5*0.5,'CData',cmap5)
end
