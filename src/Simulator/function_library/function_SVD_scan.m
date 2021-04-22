function [S_matrix]=function_SVD_scan(imageStack,x_val,y_val,...
    data_window_lat,interpol_window_lat, data_window_weights_option, G_std,...
    G_PSF, x_foc,y_foc,test_points_per_pixel)

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


S_matrix=zeros(data_window_lat^2,length(y_val),length(x_val));

numberOfImages=size(imageStack,3);

%%{
% parfor y_ind=1:length(y_val);
for y_ind=1:length(y_val);
    y=y_val(y_ind);
    S=zeros(data_window_lat^2,length(x_val));
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
        
        Gaussian_pad_=reshape(Gaussian_pad,data_window,data_window);
        Gaussian_pad_=Gaussian_pad_(Py_min:Py_max,Px_min:Px_max);
        
        if strcmpi(data_window_weights_option,'empty')
            data_window_weights=[];
        else
            data_window_weights=Gaussian_pad_(:);
        end
        [s]=function_SVD (im_, data_window_weights);
        if length(s(:)) < data_window_lat^2 
            s=[s(:);zeros(data_window_lat^2 -length(s(:)),1)];
        end
        S(:,x_ind)=s;  
 
        
    end
    S_matrix(:,y_ind,:)=S;
end
S_matrix=reshape(S_matrix,size(S_matrix,1),[]);

end