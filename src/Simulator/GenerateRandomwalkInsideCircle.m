%%
% © Arif Ahmed Sekh, UiT The Arctic University of Norway
% A Ahmed, I-S. Opstad, AB Birgisdottir, T Myrmel, BS Ahluwalia, K Agarwal, Dilip K Prasad, 
% Learning nanoscale motion patterns of vesicles in living cells, 
% IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2020
%%
clear all;
close all;
clc;
set(0,'defaultaxesfontsize',20);
set(0,'DefaultFigureVisible','off')

%parameters;
addpath('./function_library');
%addpath('./external_function_library');
%% generate string to save data
current=(clock).';
current(end)=[];
file_str_basic=['RandomWalk_'];
underscore='_';
save_path1='.\Data2\';

%% Generator property

%lysosome size (nm)
max_size=200;
min_size=100;

%r of movement (mm)
maxv=100;%velocity
maxr=300;%Circle
%number of samples need to generate

number_of_sample=2000;

flag_emission_same=false;
for s = 1700:number_of_sample
    parameters;
    sample_size=randi([min_size max_size]);
    Sample.Set(set_no).Region.params=1e-9*[sample_size];
    mkdir(strcat(save_path1,num2str(s)));
    save_path=strcat(save_path1,num2str(s),'\');
    cv=randi([0 maxv]);
    cr=randi([0 maxr]);
    Sample.Set(set_no).Dynamics=Dynamics(14);
    Sample.Set(set_no).Dynamics.params=[cv cr];
    Params.cr=cr;
    Params.cv=cv;
    Params.sample_size=sample_size;
    sim_14;
    outputFileName=[save_path file_str_basic underscore 'Image' '.tif'];
    add_poisson_noise_image3(outputFileName,250,80,s)
end
