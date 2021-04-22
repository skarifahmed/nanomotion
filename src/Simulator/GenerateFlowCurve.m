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
file_str_basic=['FlowCurve_'];
underscore='_';
save_path1='.\Data\';

%% Generator property

%lysosome size (nm)
max_size=200;
min_size=100;

%max velocity in x and y (nm)
minx=0;
maxx=10;
miny=0
maxy=10;

mina=0;
maxa=10;
minb=0;
maxb=10;
%number of samples need to generate

number_of_sample=2000;

flag_emission_same=false;
for s = 1:number_of_sample
    parameters;
    sample_size=randi([min_size max_size]);
    vx=randi([minx maxx]);
    vy=randi([miny maxy]);
    
    a=randi([mina maxa]);
    b=randi([minb maxb]);
    Sample.Set(set_no).Region.params=1e-9*[sample_size];
    mkdir(strcat(save_path1,num2str(s)));
    save_path=strcat(save_path1,num2str(s),'\');
   
    Sample.Set(set_no).Dynamics=Dynamics(9);
    Sample.Set(set_no).Dynamics.params=[vx vy 0 a b];
    Params.vx=vx;
    Params.vy=vy;
    Params.a=a;
    Params.b=b;
    Params.sample_size=sample_size;
    
    p1=randi([70 90]);
    p2=randi([200 300]);
    Params.p1=p1;
    Params.p2=p2;
    
    sim_11;
    outputFileName=[save_path file_str_basic underscore 'Image' '.tif'];
    add_poisson_noise_image1(outputFileName,p2,p1,s)
end
