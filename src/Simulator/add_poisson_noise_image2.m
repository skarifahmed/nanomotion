function add_poisson_noise_image1(im_name,max_sig,offset,fNum,varargin)
if nargin==5
    imageStack=varargin{1};
else
    ext=fliplr(strtok(fliplr(im_name),'.'));
    if strcmpi(ext,'mat')
        load(im_name);
    elseif strcmpi(ext(1:3),'tif')
        info = imfinfo(im_name);
        imageStack = [];
        numberOfImages = length(info);
        for k = 1:numberOfImages
            currentImage = imread(im_name, k, 'Info', info);
            currentImage=im2double(currentImage);
            imageStack(:,:,k) = currentImage;
        end
%         max_sig=max(imageStack(:));
    end
end
[~,rem]=strtok(fliplr(im_name),'.');
im_name=['data2/final/' num2str(fNum)];
imageStack=double(imageStack);
imageStack=uint16((max_sig-offset)*imageStack/max(imageStack(:))+offset);
imageStack=imnoise(imageStack,'poisson');
save(im_name)
im_name=[im_name '.tif'];
%%
for K = 1:size(imageStack,3)
    if K==1
        imwrite(imageStack(:, :, K), im_name, 'WriteMode', 'overwrite', 'Compression','none');
    else
        imwrite(imageStack(:, :, K), im_name, 'WriteMode', 'append', 'Compression','none');
    end
end