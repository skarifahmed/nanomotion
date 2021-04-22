%% Read me
% This file is used to setup the parameters of the optical system, sample,
% and the algorithm(s).
% This file should be used as initialization/parametrization file only and
% computations of any form should be avoided here.

% Intellectual property
% This code is created by Krishna Agarwal (BioSym, SMART, MIT) on 11th May
% 2015. The copyright of the code is held by Krishna Agarwal and SMART. The
% code can be used freely for non-commercial research purposes. Kindly cite
% the following paper if you publish the results generated using this code and
% acknowledge the copyright holders.
% !!!!!

% Essential information
% All length and time measurements are in SI units
% z-axis is the longitudinal direction

% At present it is assumed that the emitters are in the focal region only.
% the focus of the CCD is correctly set and other optical aberrations,
% misalignements, etc. are absent. Thus, just a 2 D expresssion of the PSF
% should suffice. If the things need to be changed, also change the related
% functions and scripts accordingly.

%% Measurement system parameters - Microscope

Sensor.pixel_sensor = 80e-9; % dimension of each pixel in meter
Sensor.N_x_sensor = 25;%25; % number of pixels in the x-direction
Sensor.N_y_sensor = 25;%25; % number of pixels in the y-direction
Sensor.timestep=5e-3; % time step of acquisitions
Sensor.no_of_frames=40; % total number of frames
Sensor.acquisition_time_offset=0 ; % time passed in seconds between the beginning of the dynamics of the sample and the beginning of acquisition 
Sensor.magnification = 1; % Optical magnification of the system
Sensor.wavelength = 666e-3; % wavelength of observation in micrometer.
Sensor.NA = 1.42; %numerical aperture of the system
Sensor.n_obj_n_ccd=[1.51 1]; % refractive indices of the focal and ccd region. For oil immersion, first value is the refractive index of oil. The second value, n_ccd is most often 1
Sensor.f_obj_f_ccd=[1 Sensor.magnification/1.51]*1e6;% focal lengths of the objective and ccd in microns. Generally, f_ccd/f_obj=magnification
Sensor.PSF=[];% string for special PSFs such as paraxial PSF, 2-D PSF, cofocal PSF etc.

% Sensor.Noise.Name='Gaussian';% mean of Gaussian is zero and variance
%       depends on desired SNR...If poisson of very large observations is
%       considered, then sigma=sqrt(photons per step) may be used

Sensor.Noise.Name='poisson';% lambda of poisson is sqrt(photons per step)

%% Sample information
Sample.timestep = 5e-4; %Sample.timestep: time step in seconds % this may 
%      be different from the time step of acquisition - depends on
%      fluorescent dye and quenching agents
Sample.ObservationTime=Sensor.timestep*Sensor.no_of_frames;1; % total observation time in seconds
Sample.xyz_span=[5 5 1e-3]*1e-6;% 1-D matrix containing 3 elements, the spans of sample domain - in meters

load('Libraries')
set_no=0;
%{
% DO NOT USE THIS....
% 2 emitters separated by 30 nm (at y=-15 nm)
set_no=set_no+1;
Sample.Set(set_no).NoOfEmitters=4;
Sample.Set(set_no).Region=Region(1);
Sample.Set(set_no).Region.params=1e-9*[0 0 0 15*sqrt(2)];
Sample.Set(set_no).Region.DOF=1;
Sample.Set(set_no).Pattern=Pattern(4);
Sample.Set(set_no).Pattern.min_dist=30e-9;0.9999*2*pi*Sample.Set(set_no).Region.params(4)/Sample.Set(set_no).NoOfEmitters;
Sample.Set(set_no).Blinking=Blinking(3);
Sample.Set(set_no).Blinking.default=[1e5 0.3 1e-3 0.1];
Sample.Set(set_no).Dynamics=Dynamics(1);
%}
%{
set_no=set_no+1;
Sample.Set(set_no).NoOfEmitters=60;
Sample.Set(set_no).Region=Region(1);
Sample.Set(set_no).Region.params=1e-9*[0 0 0 500];
Sample.Set(set_no).Pattern=Pattern(4);
Sample.Set(set_no).Pattern.min_dist=0.9999*2*pi*Sample.Set(set_no).Region.params(4)/Sample.Set(set_no).NoOfEmitters;
Sample.Set(set_no).Blinking=Blinking(3);
Sample.Set(set_no).Dynamics=Dynamics(1);
%}
%{
set_no=set_no+1;
Sample.Set(set_no).NoOfEmitters=180;
Sample.Set(set_no).Region=Region(1);
Sample.Set(set_no).Region.params=1e-9*[-550 0 0 500];
Sample.Set(set_no).Pattern=Pattern(4);
Sample.Set(set_no).Pattern.min_dist=0.9999*2*pi*Sample.Set(set_no).Region.params(4)/Sample.Set(set_no).NoOfEmitters;
Sample.Set(set_no).Blinking=Blinking(1);
Sample.Set(set_no).Dynamics=Dynamics(1);
%}
%{
set_no=set_no+1;
Sample.Set(set_no).NoOfEmitters=180;
Sample.Set(set_no).Region=Region(1);
Sample.Set(set_no).Region.params=1e-9*[0 550 0 500];
Sample.Set(set_no).Pattern=Pattern(4);
Sample.Set(set_no).Pattern.min_dist=0.9999*2*pi*Sample.Set(set_no).Region.params(4)/Sample.Set(set_no).NoOfEmitters;
Sample.Set(set_no).Blinking=Blinking(1);
Sample.Set(set_no).Dynamics=Dynamics(1);
%}
%{
set_no=set_no+1;
Sample.Set(set_no).NoOfEmitters=180;
Sample.Set(set_no).Region=Region(1);
Sample.Set(set_no).Region.params=1e-9*[0 -550 0 500];
Sample.Set(set_no).Pattern=Pattern(4);
Sample.Set(set_no).Pattern.min_dist=0.9999*2*pi*Sample.Set(set_no).Region.params(4)/Sample.Set(set_no).NoOfEmitters;
Sample.Set(set_no).Blinking=Blinking(1);
Sample.Set(set_no).Dynamics=Dynamics(1);
%}
%{
set_no=set_no+1;
Sample.Set(set_no).NoOfEmitters=180;
Sample.Set(set_no).Region=Region(1);
Sample.Set(set_no).Region.params=1e-9*[550 0 0 500];
Sample.Set(set_no).Pattern=Pattern(4);
Sample.Set(set_no).Pattern.min_dist=0.9999*2*pi*Sample.Set(set_no).Region.params(4)/Sample.Set(set_no).NoOfEmitters;
Sample.Set(set_no).Blinking=Blinking(1);
Sample.Set(set_no).Dynamics=Dynamics(1);
%}
%{
set_no=set_no+1;
Sample.Set(set_no).NoOfEmitters=180;
Sample.Set(set_no).Region=Region(1);
Sample.Set(set_no).Region.params=1e-9*[550 0 0 500];
Sample.Set(set_no).Pattern=Pattern(4);
Sample.Set(set_no).Pattern.min_dist=0.9999*2*pi*Sample.Set(set_no).Region.params(4)/Sample.Set(set_no).NoOfEmitters;
Sample.Set(set_no).Blinking=Blinking(1);
Sample.Set(set_no).Dynamics=Dynamics(1);
%}
%{
%Fork Arm1
set_no=set_no+1;
Sample.Set(set_no).NoOfEmitters=200;
Sample.Set(set_no).Region=Region(13);
Sample.Set(set_no).Region.params=1e-9*[0 0 500*cosd(75) 500*sind(75)];
Sample.Set(set_no).Pattern=Pattern(2);
Sample.Set(set_no).Blinking=Blinking(3);
Sample.Set(set_no).Blinking.params=[1e5 2e-3 1 0.1];
Sample.Set(set_no).Dynamics=Dynamics(1);
Sample.Set(set_no).Dynamics.params=Sample.Set(set_no).Dynamics.default;
%}
%{
%Fork Arm2
set_no=set_no+1;
Sample.Set(set_no).NoOfEmitters=200;
Sample.Set(set_no).Region=Region(13);
Sample.Set(set_no).Region.params=1e-9*[0 0 500*cosd(90+15) 500*sind(90+15)];
Sample.Set(set_no).Pattern=Pattern(2);
Sample.Set(set_no).Blinking=Blinking(3);
Sample.Set(set_no).Blinking.params=[1e5 2e-3 1 0.1];
Sample.Set(set_no).Dynamics=Dynamics(1);
Sample.Set(set_no).Dynamics.params=Sample.Set(set_no).Dynamics.default;
%}
%{
%Fork Stem;
set_no=set_no+1;
Sample.Set(set_no).NoOfEmitters=200;
Sample.Set(set_no).Region=Region(13);
Sample.Set(set_no).Region.params=1e-9*[0 0 500*cosd(-90) 500*sind(-90)];
Sample.Set(set_no).Pattern=Pattern(2);
Sample.Set(set_no).Blinking=Blinking(3);
Sample.Set(set_no).Blinking.params=[1e9 0.2 1 0.00001];
% Sample.Set(set_no).Dynamics=Dynamics(1);
% Sample.Set(set_no).Dynamics.params=Sample.Set(set_no).Dynamics.default;
Sample.Set(set_no).Dynamics=Dynamics(2);
Sample.Set(set_no).Dynamics.params=[1e-10 0.2e-5 0 0 0 0];
%}

%{
% Arif Working
set_no=set_no+1;
Sample.Set(set_no).NoOfEmitters=500;
Sample.Set(set_no).Region=Region(13);
Sample.Set(set_no).Pattern=Pattern(2);
Sample.Set(set_no).Blinking=Blinking(3);
Sample.Set(set_no).Dynamics=Dynamics(1);
%Sample.Set(set_no).Dynamics.params=Sample.Set(set_no).Dynamics.default;
Sample.Set(set_no).Dynamics=Dynamics(2);
Sample.Set(set_no).Dynamics.params=[1e-10 0.2e-5 0 0 0 0];
%}

%Arif new
set_no=set_no+1;
Sample.Set(set_no).NoOfEmitters=1000;
Sample.Set(set_no).Region=Region(14);
%%Size
Sample.Set(set_no).Region.params=1e-9*[100];
Sample.Set(set_no).Pattern=Pattern(2);
Sample.Set(set_no).Blinking=Blinking(3);

%stationary Lysosome
%Sample.Set(set_no).Dynamics=Dynamics(1);
%Sample.Set(set_no).Dynamics.params=Sample.Set(set_no).Dynamics.default;

%Moving Lysosome
%Sample.Set(set_no).Dynamics=Dynamics(2);
%Sample.Set(set_no).Dynamics.params=[1e-10 0.2e-5 0 0 0 0];

%Random Walk Lysosome
%Sample.Set(set_no).Dynamics=Dynamics(6);
%Sample.Set(set_no).Dynamics.params=[1 1 0];

%Random Walk Lysosome in a Circle
%Sample.Set(set_no).Dynamics=Dynamics(10);
%Sample.Set(set_no).Dynamics.params=[50];


%Random Walk Lysosome in a Circular Moton
%Sample.Set(set_no).Dynamics=Dynamics(11);
%Sample.Set(set_no).Dynamics.params=[0.2];

%Random Walk Lysosome in a Circular Moton r Range
%Sample.Set(set_no).Dynamics=Dynamics(12);
%Sample.Set(set_no).Dynamics.params=[10 50];

%Random Walk Lysosome inside circle
%Sample.Set(set_no).Dynamics=Dynamics(13);
%Sample.Set(set_no).Dynamics.params=[5 300];

%Random Walk Lysosome Small region
Sample.Set(set_no).Dynamics=Dynamics(11);
Sample.Set(set_no).Dynamics.params=[50 200];

%Random Walk Lysosome Curve
%Sample.Set(set_no).Dynamics=Dynamics(9);
%Sample.Set(set_no).Dynamics.params=[6 10 0 6 4];

%%Random Flow Lysosome
%Sample.Set(set_no).Dynamics=Dynamics(6);
%Sample.Set(set_no).Dynamics.params=[1e-10 0.2e-5 0 0 0 0];


%%Flow and Stationary Lysosome
%Sample.Set(set_no).Dynamics=Dynamics(7);
%Sample.Set(set_no).Dynamics.params=[1e-10 0.2e-5 0 0 0 0];

%%Flow and Random Lysosome
%Sample.Set(set_no).Dynamics=Dynamics(8);
%Sample.Set(set_no).Dynamics.params=[1e-10 0.2e-5 0 0 0 0 3 4 0];

%Sample.Set(set_no).Dynamics=Dynamics(2);
%Sample.Set(set_no).Dynamics.params=[1e-10 0.2e-5 0 0 0 0];


%{
% Five lines
for position=(-100:100:100)
set_no=set_no+1;
Sample.Set(set_no).NoOfEmitters=11;
Sample.Set(set_no).Region=Region(13);
Sample.Set(set_no).Region.params=[position -250 position 250]*1e-9;
Sample.Set(set_no).Pattern=Pattern(4);
Sample.Set(set_no).Blinking=Blinking(3);
Sample.Set(set_no).Dynamics=Dynamics(1);
Sample.Set(set_no).Dynamics.params=Sample.Set(set_no).Dynamics.default;
end
%}
%{
set_no=set_no+1;
Sample.Set(set_no).NoOfEmitters=1;
Sample.Set(set_no).Region=Region(1);
Sample.Set(set_no).Region.params=1e-6*[-0.2 0 0 0.05]*3;
Sample.Set(set_no).Region.DOF=2;
Sample.Set(set_no).Pattern=Pattern(3);
Sample.Set(set_no).Pattern.params=true;
Sample.Set(set_no).Blinking=Blinking(1);
Sample.Set(set_no).Dynamics=Dynamics(2);
Sample.Set(set_no).Dynamics.params=1e-6 * [0 1 0 0 0 0]*3;
%}

%% spokes pattern
%{
for set_no=1:5%36
Sample.Set(set_no).NoOfEmitters=201;
Sample.Set(set_no).Region=Region(13);
Sample.Set(set_no).Region.params=[50*cos(set_no*pi/18) 50*sin(set_no*pi/18) 700*cos(set_no*pi/18) 700*sin(set_no*pi/18)]*1e-9;
Sample.Set(set_no).Pattern=Pattern(4);
Sample.Set(set_no).Blinking=Blinking(3);
Sample.Set(set_no).Blinking.params=Blinking(3).default;
Sample.Set(set_no).Blinking.params(2)=0.0003;
Sample.Set(set_no).Blinking.params(2)=1;
Sample.Set(set_no).Dynamics=Dynamics(1);
Sample.Set(set_no).Dynamics.params=Sample.Set(set_no).Dynamics.default;
end
%}
%{
% Lines in 1 micron by 1 micron region
% lpmicron=15;
x=linspace(-500,500,lpmicron+1);
p_blink=0.0003;tau_blink=1;
for set_no=1:(lpmicron+1)
Sample.Set(set_no).NoOfEmitters=201;
Sample.Set(set_no).Region=Region(13);
Sample.Set(set_no).Region.params=[x(set_no) -500 x(set_no) +500]*1e-9;
Sample.Set(set_no).Pattern=Pattern(4);
Sample.Set(set_no).Blinking=Blinking(3);
Sample.Set(set_no).Blinking.params=Blinking(3).default;
Sample.Set(set_no).Blinking.params(2)=p_blink;
Sample.Set(set_no).Blinking.params(2)=tau_blink;
Sample.Set(set_no).Dynamics=Dynamics(1);
Sample.Set(set_no).Dynamics.params=Sample.Set(set_no).Dynamics.default;
end
%}
%{
% grid of 100 emitters placed in 10 x 10 configuration
% lpmicron=15;
N=10;
x=N*1000*linspace(-0.5*(1-1/N),0.5*(1-1/N),N)/lpmicron;
p_blink=0.0003;tau_blink=1;
for set_no=1:N
Sample.Set(set_no).NoOfEmitters=10;
Sample.Set(set_no).Region=Region(13);
Sample.Set(set_no).Region.params=[x(set_no) x(1) x(set_no) +x(N)]*1e-9;
Sample.Set(set_no).Pattern=Pattern(4);
Sample.Set(set_no).Blinking=Blinking(3);
Sample.Set(set_no).Blinking.params=Blinking(3).default;
Sample.Set(set_no).Blinking.params(2)=p_blink;
Sample.Set(set_no).Blinking.params(2)=tau_blink;
Sample.Set(set_no).Dynamics=Dynamics(1);
Sample.Set(set_no).Dynamics.params=Sample.Set(set_no).Dynamics.default;
end
%}
%{
% Pair of lines separated by given number
% lpmicron=15;
x=1000*[-1 1]*(0.5/lpmicron); % linspace(-500,500,lpmicron+1);
p_blink=0.0003;tau_blink=1;
for set_no=1:2
Sample.Set(set_no).NoOfEmitters=51;
Sample.Set(set_no).Region=Region(13);
Sample.Set(set_no).Region.params=[x(set_no) -250 x(set_no) +250]*1e-9;
Sample.Set(set_no).Pattern=Pattern(4);
Sample.Set(set_no).Blinking=Blinking(3);
Sample.Set(set_no).Blinking.params=Blinking(3).default;
Sample.Set(set_no).Blinking.params(2)=p_blink;
Sample.Set(set_no).Blinking.params(2)=tau_blink;
Sample.Set(set_no).Dynamics=Dynamics(1);
Sample.Set(set_no).Dynamics.params=Sample.Set(set_no).Dynamics.default;
end
%}




