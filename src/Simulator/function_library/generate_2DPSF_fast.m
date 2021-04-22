%% Read me
% This function generates the 2-D PSF of the system very fast.
% It computes the PSF for minimum required combinations, and then pads and 
% shifts these combinations to get the PSF fast.
% It may be made to use a PSF computation algorithm of desire, for example,
% generate_PSF, generate_PSF_Airy, generate_PSF_gaussian, etc.

% Intellectual property
% This code is created by Krishna Agarwal (BioSym, SMART, MIT) on 5th Dec
% 2014. The copyright of the code is held by Krishna Agarwal and SMART. The
% code can be used freely for non-commercial research purposes. Kindly cite
% the paper below if you publish the results generated using this code and
% acknowledge the copyright holders.
%!!!!!!!!!!

% Inputs:
% x_sensor: x-coordinates of the sensor pixels, N_x_pixels x N_y_pixels
% y_sensor: y-coordinates of the sensor pixels, N_x_pixels x N_y_pixels
% x_foc: x-coordinates of the focal domain, N_x_focalregions x N_y_focalregions
% y_foc: y-coordinates of the focal domain, N_x_focalregions x N_y_focalregions
% Sensor: structure containing at least the follwoing fields (defined in
%   parameters.m)
%   pixel_sensor: dimension of each pixel in micrometer
%   magnification: Optical magnification of the system
%   NA: numerical aperture of the system
%   n_obj_n_ccd: refractive indices of the focal and ccd region. For oil immersion, first value is the refractive index of oil. The second value, n_ccd is most often 1
%   wavelength: wavelength of observation in micrometer

% Outputs:
% I_map: intensity PSF from focus [0 0 0] to the input ccd points
% G: dyadic PSF from focus [0 0 0] to the input ccd points
% I_out: mapping from the input sample points to the input CCD points

function [I_map,G,I_out]=generate_2DPSF_fast(x_ccd,y_ccd,x_sam,y_sam,Sensor,num_noise_level)

%step 1a, how many min. pixels window is enough to represent PSF (min_win)
%step 1b, calculate the computable pixels
%step 2a, how many combinations of PSF need to be computed (min_comb)
%step 2b, calculate these PSFs for the min. pixel window (PSF_min_win_min_comb)
%step 3a, what is the largest ccd-foc distance encountered (max_win)
%step 3b, pad the PSF_min_win_min_comb upto the max_win (PSF_max_win_min_comb)
%step 4, shift the PSF_max_win_min_comb for each focal pixel to get the
%PSF_matrix


%% step1a
x_test=0:100*Sensor.pixel_sensor;
[phi,th,r]=cart2sph(x_test,zeros(size(x_test)),zeros(size(x_test)));
vec_r_ccd_temp=[r(:),th(:),phi(:)]; clear r th phi;
if isempty(Sensor.PSF)
    [~,~,PSF_est]=generate_PSF(Sensor,[0 0 0],vec_r_ccd_temp);
else
    eval(['[~,G,PSF_est]=generate_PSF_' Sensor.PSF '(Sensor,[0 0 0],vec_r_ccd_temp);']);
end
min_win=sum(PSF_est>num_noise_level*max(PSF_est))+1;
x_min_win=Sensor.pixel_sensor*(-min_win:min_win);
y_min_win=Sensor.pixel_sensor*(-min_win:min_win);
% [x_min_win,y_min_win]=meshgrid(x_min_win,y_min_win);
% x_min_win=x_min_win(:);y_min_win=y_min_win(:);


%% step 1b
[x_sensor_grid,y_sensor_grid]=meshgrid(x_ccd,y_ccd);
[x_foc_grid,y_foc_grid]=meshgrid(x_sam,y_sam);

x_sensor=x_sensor_grid(:);
y_sensor=y_sensor_grid(:);
x_foc=x_foc_grid(:);
y_foc=y_foc_grid(:);

%% step2a

sensor_pixel_in_foc=Sensor.pixel_sensor/Sensor.magnification;
% 
% [xx_sensor_ind,xx_foc_ind]=meshgrid(x_sensor/Sensor.pixel_sensor,x_foc/sensor_pixel_in_foc);
% [yy_sensor_ind,yy_foc_ind]=meshgrid(y_sensor/Sensor.pixel_sensor,y_foc/sensor_pixel_in_foc);
% 
% [tt,dd]=cart2pol(xx_sensor_ind-xx_foc_ind,yy_sensor_ind-yy_foc_ind);
% dd=(1e-6)*round(1e6 *dd);
% figure;imagesc(dd)
% figure;imagesc(tt)
% [ddd,iddd]=min(dd,[],2);
% idddd=sub2ind(size(tt),(1:size(dd,1)).',iddd);
% ttt=tt(idddd);
% [C1,ia1,ic1]=unique([ddd(:),ttt(:)],'rows');

x_foc_rel_pixel=(1e-6)*round(1e6 * rem(x_foc,sensor_pixel_in_foc));
y_foc_rel_pixel=(1e-6)*round(1e6 * rem(y_foc,sensor_pixel_in_foc));

[C,ia,ic]=unique([x_foc_rel_pixel(:) y_foc_rel_pixel(:)],'rows');

%% step 2b
[xx,yy]=meshgrid(x_min_win,y_min_win);
[phi,th,r]=cart2sph(xx,yy,zeros(size(xx)));clear xx yy;
vec_r_ccd=[r(:),th(:),phi(:)]; clear r th phi;
[phi,th,r]=cart2sph(C(:,1),C(:,2),zeros(size(C(:,1))));
vec_r_foc=[r(:),th(:),phi(:)]; clear r th phi;
if isempty(Sensor.PSF)
    [I_map,G,PSF_min_win_min_comb]=generate_PSF(Sensor,vec_r_foc,vec_r_ccd);
else
    eval(['[I_map,G,PSF_min_win_min_comb]=generate_PSF_' Sensor.PSF '(G,vec_r_foc,vec_r_ccd);']);
end
PSF_min_win_min_comb=reshape(PSF_min_win_min_comb,2*min_win+1,2*min_win+1,[]);
%% step 3a
[xx_sensor_ind,xx_foc_ind]=meshgrid(x_sensor/Sensor.pixel_sensor,x_foc/sensor_pixel_in_foc);
[yy_sensor_ind,yy_foc_ind]=meshgrid(x_sensor/Sensor.pixel_sensor,x_foc/sensor_pixel_in_foc);
max_win=max(min_win,ceil(max(max(max(abs(xx_sensor_ind-xx_foc_ind))),max(max(abs(yy_sensor_ind-yy_foc_ind))))));

%% step 3b
PSF_max_win_min_comb=zeros(2*max_win+1,2*max_win+1,length(Cx));
PSF_max_win_min_comb(max_win+1+(-min_win:min_win),max_win+1+(-min_win:min_win),:)=PSF_min_win_min_comb;
PSF_max_win_min_comb=reshape(PSF_max_win_min_comb,(2*max_win+1).^2,[]);
%% step 4
x_ind_ccd=(max_win+1)+(xx_sensor_ind-xx_foc_ind);
y_ind_ccd=(max_win+1)+(yy_sensor_ind-yy_foc_ind);

inf_ccd=sub2ind([(2*max_win+1) (2*max_win+1)],x_ind_ccd,y_ind_ccd);
I_out=PSF_max_win_min_comb(x_ind_ccd,y_ind_ccd,ic);

return