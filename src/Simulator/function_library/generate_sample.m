%% Read me
% This function generates a distribution of emitters and their blinking
% statsitics using the properties of emitters and their blinking dynamics.

% Intellectual property
% This code is created by Krishna Agarwal (BioSym, SMART, MIT) on 26th Sept
% 2014. The copyright of the code is held by Krishna Agarwal and SMART. The
% code can be used freely for non-commercial research purposes. Kindly cite
% the paper below if you publish the results generated using this code and
% acknowledge the copyright holders.
%!!!!!!!!!!

% Inputs:
% Sample: a structure containing the following elements:
%     Sample.N_emitters : number of emitters in the sample
%     Sample.r_min_emitters : minimum radius of emitter in nanometers
%     Sample.r_max_emitters : maximum radius of emitter in nanometers
%     Sample.dist_min_edges_emitters : minimum distance between the edges of two closest emitters in nanometers
%     Sample.dist_max_edges_emitters : maximum distance between the edges of two closest emitters in nanometers
%     Sample.emitter_strength_min : the min strength of the emitter, it should be the dipole moment of the emitter. Ohm/meter (check!)
%     Sample.emitter_strength_max : the max strength of the emitter, it should be the dipole moment of the emitter. Ohm/meter (check!)
%     Sample.time_series : time series of sample's emission dynamics in milliseconds.
%     Sample.PDF = [Tau_ON_mean P_ON]
%       Tau_ON_mean (millisecond) - average mean time determined using Possion distribution of the intensity
%       P_ON: probability of an emitter being ON at a particulr time.
% distribution_map: probability distribution map of the presence of
%                   emitter.
% x_foc,y_foc,z_foc : 1D vectors containing the x,y,z-coordinates of the
%                       focal region's computation grid.
% distribution_map_name: Name of the distrbution (important for Npointxxx
%       types of distribution

% Outputs:
% Stack: the actual emission blinking dynamics (eps_k*s_k(t)
%       of Dertinger 2009)
%       4-D, N_x_foc x N_y_foc x N_z_foc x N_t_foc
% Emission: the emission strengths (eps_k of Dertinger 2009)
%       3-D, N_x_foc x N_y_foc x N_z_foc
% time_series: 1-D matrix containing the time coordinates of the
%       blinking timeline
%       N_t_foc x 1
% N_t_foc: number of time samples in time_series


function [Stack,Emission,time_series,N_t_foc]=generate_sample(Sample,distribution_map,x_foc,y_foc,z_foc,distribution_map_name)
%% inputs initialization
N_emitters=Sample.N_emitters;
r_min_emitters=Sample.r_min_emitters;
r_max_emitters=Sample.r_max_emitters;
dist_min_edges_emitters=Sample.dist_min_edges_emitters;
if isempty(dist_min_edges_emitters)
    dist_min_edges_emitters=0;
end
dist_max_edges_emitters=Sample.dist_max_edges_emitters;
if isempty(dist_max_edges_emitters)
    dist_max_edges_emitters=sqrt((x_foc(end)-x_foc(1)).^2 + (y_foc(end)-y_foc(1)).^2 + (z_foc(end)-z_foc(1)).^2);
end
emitter_strength_min=Sample.emitter_strength_min;
if isempty(emitter_strength_min)
    emitter_strength_min=min(emitter_strength_max,1);
end
emitter_strength_max=Sample.emitter_strength_max;
if isempty(emitter_strength_max)
    emitter_strength_max=emitter_strength_min;
end
time_series=Sample.time_series;
N_t_foc=length(time_series);

Tau_ON_mean=Sample.PDF(1); % millisecond - average mean time determined using Possion distribution of the intensity
Prob_ON=Sample.PDF(2);

clear Sample;
%% Outputs and temp initializations
Emission=zeros(size(distribution_map));
Stack=zeros([size(distribution_map),N_t_foc]);

radii= r_min_emitters + rand([1,N_emitters])*(r_max_emitters-r_min_emitters);
emission_strengths= emitter_strength_min + rand([1,N_emitters])*(emitter_strength_max-emitter_strength_min);

T_chunk=Tau_ON_mean/(time_series(2)-time_series(1));
%% choosing centers of emitters
[xx,yy,zz]=meshgrid(x_foc,y_foc,z_foc);
found_solution=false;
count=0;
while(~found_solution)
    count=count+1;
    if length(distribution_map_name)>6 & strcmpi(distribution_map_name(1:6),'Npoint')
        length(find(distribution_map))==2;
        center_indices=find(distribution_map);
        found_solution= true;
        N_emitters=length(center_indices);
        radii=radii(1)*ones(1,N_emitters);
        emission_strengths=emission_strengths(1)*ones(1,N_emitters);
    else
        center_indices = gendist(distribution_map(:)/sum(distribution_map(:)),1,N_emitters);
        if N_emitters==1
            ends_x=[1 1 1 1 length(xx(:)) length(xx(:)) length(xx(:)) length(xx(:))];
            ends_y=[1 1 length(yy(:)) length(yy(:)) 1 1 length(yy(:)) length(yy(:))];
            ends_z=[1 length(zz(:)) 1 length(zz(:)) 1 length(zz(:)) 1 length(zz(:))];
            dist=sqrt((xx(ends_x)-xx(center_indices)).^2 + (yy(ends_y)-yy(center_indices)).^2 + (zz(ends_z)-zz(center_indices)).^2);
            if min(dist)>radii
                found_solution= true;
            elseif count==100000
                error('Could not generate sample');
                break;
            end
        else
            [xx1,xx2]=meshgrid(xx(center_indices),xx(center_indices));
            [yy1,yy2]=meshgrid(yy(center_indices),yy(center_indices));
            [zz1,zz2]=meshgrid(zz(center_indices),zz(center_indices));
            [r1,r2]=meshgrid(radii,radii);
            dist=sqrt((xx1-xx2).^2+(yy1-yy2).^2+(zz1-zz2).^2) - (r1+r2);
            dist=dist-diag(diag(dist));
            %     max_dist=max(dist(:));
            dist(dist == 0)=Inf;
            min_dist = min(dist(:));
            max_dist = max(min(dist));
            if min_dist>dist_min_edges_emitters && max_dist<dist_max_edges_emitters
                found_solution= true;
            elseif count==100000
                error('Could not generate sample');
                break;
            end
        end
    end
end

%% Computing outputs
Stack=reshape(Stack,[],N_t_foc);

for n=1:N_emitters
    center_ind=center_indices(n);
    ttt=[];
    if length(x_foc)>1
        ttt=[ttt;(x_foc(1)-x_foc(2)).^2];
    end
    if length(y_foc)>1
        ttt=[ttt;(y_foc(1)-y_foc(2)).^2];
    end
    if length(z_foc)>1
        ttt=[ttt;(z_foc(1)-z_foc(2)).^2];
    end
    
    indices = round(((xx-xx(center_ind)).^2+(yy-yy(center_ind)).^2+(zz-zz(center_ind)).^2 - radii(n).^2)/min(ttt))<=0;
    Emission(indices)=emission_strengths(n);
    
    no_of_ONs=random('Binomial',T_chunk,Prob_ON,[1 floor(N_t_foc/T_chunk)]);
    blinking=zeros(1,N_t_foc);
    for nn=1:floor(N_t_foc/T_chunk)
        if no_of_ONs(nn)~=0
            if no_of_ONs(nn)~=T_chunk
                tt=randperm(T_chunk-no_of_ONs(nn),1);
            else
                tt=0;
            end
            blinking((tt+(nn-1)*T_chunk+1):(tt+(nn-1)*T_chunk+no_of_ONs(nn)))=1;
        end
    end
    Stack(indices,:) = emission_strengths(n)*repmat(blinking,[length(find(indices)) 1]);
end
Stack=reshape(Stack,[size(distribution_map),N_t_foc]);

end