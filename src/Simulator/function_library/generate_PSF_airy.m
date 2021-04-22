%% Read me
% This function generates the 2-D gaussian approximation of the PSF of the
% system.
% It has two options
% Option A: Input is a 2 element array, 
%           computes Input(1)*exp((x.^2 + y.^2)/(Input(2).^2))
% Option B: use the actual PSF of the system to first determine a and s for
%           a*exp((x.^2 + y.^2)/(s.^2)), Input is the Sensor structure made
%           in parameters.m
% Other simplified versions have also been coded - see functions
% generate_PSF_<xxxx>

% Intellectual property
% This code is created by Krishna Agarwal (BioSym, SMART, MIT) on 7th Nov
% 2014. The copyright of the code is held by Krishna Agarwal and SMART. The
% code can be used freely for non-commercial research purposes. Kindly cite
% the paper below if you publish the results generated using this code and
% acknowledge the copyright holders.
%!!!!!!!!!!

% Inputs:
% Input: first input defining the type of computation
%   If option A, Input is 2 element array [a s] used to compute 
%   a*exp((x.^2 + y.^2)/(s.^2))
%   If option B, Input is Sensor structure containing at least the 
%   following fields (defined in parameters.m)
%       NA: numerical aperture of the system
%       n_obj_n_ccd: refractive indices of the focal and ccd region. For oil immersion, first value is the refractive index of oil. The second value, n_ccd is most often 1
%       f_obj_f_ccd: focal lengths of the objective and ccd in microns. Generally, f_ccd/f_obj=magnification
%       wavelength: wavelength of observation in micrometer
% vec_r_foc: N_foc x 3 matrix: [r theta phi] cooridnates of the sample
%           points, r is in micrometers
% vec_r_ccd: N_ccd x 3 matrix: [r theta phi] cooridnates of the ccd's pixel
%           points, r is in micrometers

% Outputs:
% I_map: intensity PSF from focus [0 0 0] to the input ccd points
% G: [a s] parameters of the gaussuan function
% I_out: mapping from the input sample points to the input CCD points

function [I_map,G,I_out] = generate_PSF_gaussian(Input,vec_r_foc,vec_r_ccd)

if ~ isstruct(Input)
    lam = Input(1);% wavelength of observation in micrometer.
    M = Input(2);% Optical magnification of the system
    NA = Input(3);%numerical aperture of the system
    G=[lam M NA];
else
    Sensor = Input;
    M=Sensor.magnification; % Optical magnification of the system
    lam=Sensor.wavelength; % wavelength of observation in micrometer.
    NA=Sensor.NA; %numerical aperture of the system
    G=[lam M NA];
end

k=2*pi/lam;

[x_foc,y_foc,z_foc]=sph2cart(vec_r_foc(:,3),vec_r_foc(:,2),vec_r_foc(:,1));
[x_ccd,y_ccd,z_ccd]=sph2cart(vec_r_ccd(:,3),vec_r_ccd(:,2),vec_r_ccd(:,1));

[X_ccd,X_foc]=meshgrid(x_ccd,M*x_foc);
[Y_ccd,Y_foc]=meshgrid(y_ccd,M*y_foc);
x=X_ccd-X_foc;
y=Y_ccd-Y_foc;

rho=sqrt(x.^2 + y.^2);
rho_ = NA*rho/M;
I_out=(besselj(1,k*rho_)./(k*rho_)).^2;
I_out(isnan(I_out))=(0.5/gamma(2)).^2;


rho=sqrt(x_ccd.^2 + y_ccd.^2);
rho_ = NA*rho/M;
PSF=(besselj(1,k*rho_)./(k*rho_)).^2;
PSF(isnan(PSF))=(0.5/gamma(2)).^2;
I_map{1}=PSF/3;
I_map{2}=I_map{1};
I_map{3}=I_map{1};


