function [LSM_out,LSM_time_out]=function_LSM_scan(imageStack,x_val,y_val,...
    data_window_lat,interpol_window_lat, data_window_weights_option, G_std,...
    G_PSF, x_foc,y_foc,test_points_per_pixel, flip_data)

% data window pixels
data_pixels=sqrt(size(G_PSF,1));
ceil_half_data_window=ceil(sqrt(size(G_PSF,1))/2);
data_window_pixels=(ceil_half_data_window-floor(data_window_lat/2)):(ceil_half_data_window+floor(data_window_lat/2));
length_data_window_pixels=length(data_window_pixels);

% test interpolation window pixels
test_points=sqrt(size(G_PSF,2));
ceil_half_test_window=ceil(sqrt(size(G_PSF,2))/(2*test_points_per_pixel));
test_window_pixels=(test_points_per_pixel*(ceil_half_test_window-floor(interpol_window_lat/2) -1 ) +1):(test_points_per_pixel*(ceil_half_test_window+floor(interpol_window_lat/2)));
length_test_window_pixels=length(test_window_pixels);


% reducing PSF to desired data and test point selections
G_PSF=reshape(G_PSF,data_pixels,data_pixels,test_points,test_points);
G_PSF=G_PSF(data_window_pixels,data_window_pixels,test_window_pixels,test_window_pixels);
G_PSF=reshape(G_PSF,length_data_window_pixels^2,length_test_window_pixels^2);
x_foc=x_foc(test_window_pixels);
y_foc=y_foc(test_window_pixels);

% making filter pad;
floor_half_data_window=floor(sqrt(size(G_PSF,1))/2);
data_window=sqrt(size(G_PSF,1));
[x_temp,y_temp]=meshgrid(1:data_window,1:data_window);
Gaussian_pad=mvnpdf([x_temp(:) y_temp(:)],floor_half_data_window+[1 1],[G_std G_std]);
Gaussian_pad=Gaussian_pad/max(Gaussian_pad);

% used later in interleaving the test windows:
floor_half_interpol_window=floor(interpol_window_lat/2);
ceil_half_interpol_window=ceil(interpol_window_lat/2);

numberOfImages=size(imageStack,3);

%%{
for y_ind=1:length(y_val);
    y=y_val(y_ind);
    for x_ind=1:length(x_val);
        x=x_val(x_ind);
        
        y_min=max([y-floor_half_data_window,1]);
        y_max=min([y+floor_half_data_window,size(imageStack,1)]);
        x_min=max([x-floor_half_data_window,1]);
        x_max=min([x+floor_half_data_window,size(imageStack,2)]);
        im=imageStack(y_min:y_max,x_min:x_max,:);
        im_=reshape(im,[],numberOfImages);
        
        
        Py_min=max(1,(floor_half_data_window+1)-(y-1));
        Py_max=min(2*floor_half_data_window+1,(floor_half_data_window+1)+(size(imageStack,1)-y));
        Px_min=max(1,(floor_half_data_window+1)-(x-1));
        Px_max=min(2*floor_half_data_window+1,(floor_half_data_window+1)+(size(imageStack,2)-x));
        G_PSF_=reshape(G_PSF,data_window,data_window,[]);
        G_PSF_=G_PSF_(Py_min:Py_max,Px_min:Px_max,:);
        G_PSF_=reshape(G_PSF_,[],size(G_PSF_,3));
        
        Gaussian_pad_=reshape(Gaussian_pad,data_window,data_window);
        Gaussian_pad_=Gaussian_pad_(Py_min:Py_max,Px_min:Px_max);
        
        if strcmpi(data_window_weights_option,'empty')
            data_window_weights=[];
        elseif strcmpi(data_window_weights_option,'PSF')
            data_window_weights=G_PSF_(:,ceil(size(G_PSF_,2)/2));
        else
            data_window_weights=Gaussian_pad_(:);
        end
        [Pseudospectrum,U]=function_LSM (im_, G_PSF_,data_window_weights);
        
        current=reshape(Pseudospectrum,length(y_foc),length(x_foc));
%         current_time_reconstruct=reshape(U.',length(y_foc),length(x_foc),size(U,1));
         
        if flip_data
            current=flipdim(flipdim(current,1),2);
%             current_time_reconstruct=flipdim(flipdim(current_time_reconstruct,1),2);
        end
        
        for y_i=(-1*floor_half_interpol_window):floor_half_interpol_window
            for x_i=(-1*floor_half_interpol_window):floor_half_interpol_window
                yy=y_i+ceil_half_interpol_window;
                xx=x_i+ceil_half_interpol_window;
                im_no= (yy-1)*interpol_window_lat+xx;
                size_y_foc_pixel=round(length(y_foc)/interpol_window_lat);
                size_x_foc_pixel=round(length(x_foc)/interpol_window_lat);
                if (y_ind+y_i)>=1 && (y_ind+y_i)<=length(y_val) && (x_ind+x_i)>=1 && (x_ind+x_i)<=length(x_val)
                        LSM{y_ind+y_i,x_ind+x_i,yy,xx}=current(((yy-1)*size_y_foc_pixel+1):(yy*size_y_foc_pixel),((xx-1)*size_x_foc_pixel+1):(xx*size_x_foc_pixel));
                        mean_data_window{y_ind+y_i,x_ind+x_i,yy,xx}=ones(size_y_foc_pixel,size_x_foc_pixel);
%                             LSM_time{y_ind+y_i,x_ind+x_i,yy,xx}=current_time_reconstruct(((yy-1)*size_y_foc_pixel+1):(yy*size_y_foc_pixel),((xx-1)*size_x_foc_pixel+1):(xx*size_x_foc_pixel),:);
                end
            end
        end
        
    end
end
for y_ind=1:length(y_val);
    for x_ind=1:length(x_val);
        for xx=1:interpol_window_lat
            for yy=1:interpol_window_lat
                if isempty(LSM{y_ind,x_ind,yy,xx})
                    size_y_foc_pixel=round(length(y_foc)/interpol_window_lat);
                    size_x_foc_pixel=round(length(x_foc)/interpol_window_lat);
                    LSM{y_ind,x_ind,yy,xx}=zeros(size_y_foc_pixel,size_x_foc_pixel);
                    mean_data_window{y_ind,x_ind,yy,xx}=zeros(size_y_foc_pixel,size_x_foc_pixel);
%                     LSM_time{y_ind,x_ind,yy,xx}=zeros(size_y_foc_pixel,size_x_foc_pixel,size(imageStack,3));
                
                end
            end
        end
    end
end
%}



mean_data_window=cell2mat(mean_data_window);

LSM__=cell2mat(LSM);clear LSM;
LSM_=sum(sum(LSM__,4),3)./sum(sum(mean_data_window,4),3);clear LSM__
LSM_out=abs(1./LSM_);clear LSM_


LSM_time_out=[];
%{

for dim1=1:size(LSM_time,1)
    Row=[];
    for dim2=1:size(LSM_time,2)
        Cell=zeros([size(LSM_time{1,1,1,1})]);
        for dim3=1:size(LSM_time,3)
            for dim4=1:size(LSM_time,4)
                Cell=Cell+squeeze(cell2mat(LSM_time(dim1,dim2,dim3,dim4)));
            end
        end
        Row=[Row Cell];
    end
    LSM_time_out=[LSM_time_out;Row];
end
clear LSM_time;
%}

end