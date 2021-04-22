%% Read me
% This function generates the grid for the sensor. The centers of the
% pixels are given as the output.

% Intellectual property
% This code is created by Krishna Agarwal (BioSym, SMART, MIT) on 25th Sept
% 2014. The copyright of the code is held by Krishna Agarwal and SMART. The
% code can be used freely for non-commercial research purposes. Kindly cite
% the following paper if you publish the results generated using this code and
% acknowledge the copyright holders.
% !!!!!

% Inputs:
% Sensor.pixel_sensor : dimension of each pixel in micrometer
% Sensor.N_x_sensor : number of pixels in the x-direction
% Sensor.N_y_sensor : number of pixels in the y-direction

% Outputs:
% x_sensor,y_sensor,z_sensor : 1D vectors containing the x,y,z-coordinates
%                               of the centers of the pixels of the sensor.

function [ x_sensor,y_sensor,z_sensor ] = generate_sensor( Sensor )
z_sensor=0;
if Sensor.N_x_sensor==1
    x_sensor=0;
else
    x_sensor=linspace(-1,1,Sensor.N_x_sensor);
    x_sensor=Sensor.pixel_sensor*x_sensor/(x_sensor(2)-x_sensor(1));
end
if Sensor.N_y_sensor==1
    y_sensor=0;
else
    y_sensor=linspace(-1,1,Sensor.N_y_sensor);
    y_sensor=Sensor.pixel_sensor*y_sensor/(y_sensor(2)-y_sensor(1));
end
end

