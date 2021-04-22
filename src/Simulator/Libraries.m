%% This file contains libraries of four things:
%{
1. Region in which the emitters are restricted or inititally present
2. Distribution pattern - random or uniform distribution
3. Particle motion - stationary, diffuse, flow, or flow diffuse
4. Blinking statistics 
%}
%%
%%
%% Region library for each Distribution
%{
Region Structure contains the following fields:
    name: name of the shape
    no_of_params: number of parameters used to describe the shape
    explain_params: the explanation and sequence of parameters
    units:params: units of the parameters
    default: default values
    DOF: 1 (edge in x-y plane), 2 (surface in x-y plane), 3 (volume)
%}
%% Sphere/Circle type (x-x0).^2 + (y-y0).^2 + (z-z0).^2 <= r.^2;
Region(1).name='Sphere'; % shape
Region(1).no_of_params=4; % no. of parameters
Region(1).explain_params = {'x0','y0','z0','r'}; % explanation of parameters
Region(1).units_params={'m','m','m','m'};
Region(1).default=1e-6*[0 0 0 0.1];
Region(1).params=[];
Region(1).DOF=1; % allowed 1,2,3
%% Ellipsoid/Ellipse type ((x-x0)/a).^2 + ((y-y0)/b).^2 + ((z-z0)/c).^2 <= 1;
Region(2).name='Ellipsoid'; % shape
Region(2).no_of_params=6; % no. of parameters
Region(2).explain_params = {'x0','y0','z0','a','b','c'}; % explanation of parameters
Region(2).units_params={'m','m','m','m','m','m'};
Region(2).default=1e-6*[0 0 0 0.1 0.1 0.1];
Region(2).params=[];
Region(2).DOF=1; % allowed 1,2,3
%% Cube/Square type -a/2 <= (x-x0), (y-y0), (z-z0) <= a/2
Region(3).name='Cube'; % shape
Region(3).no_of_params=4; % no. of parameters
Region(3).explain_params = {'x0','y0','z0','a'}; % explanation of parameters
Region(3).units_params={'m','m','m','m'};
Region(3).default=1e-6*[0 0 0 0.1];
Region(3).params=[];
Region(3).DOF=1; % allowed 1,2,3
%% Cuboid/Rectangle type -a/2 <= (x-x0) <= a/2, -b/2 <= (y-y0) <= b/2, -c/2<= (z-z0) <= c/2
Region(4).name='Cuboid'; % shape
Region(4).no_of_params=6; % no. of parameters
Region(4).explain_params = {'x0','y0','z0','a','b','c'}; % explanation of parameters
Region(4).units_params={'m','m','m','m','m','m'};
Region(4).default=1e-6*[0 0 0 0.1 0.1 0.1];
Region(4).params=[];
Region(4).DOF=1; % allowed 1,2,3
%% Sphere/Circle Shell Concentric type: r_in.^2 <= (x-x0).^2 + (y-y0).^2 + (z-z0).^2 <= r_out.^2;
Region(5).name='SphereShellConcentric'; % shape
Region(5).no_of_params=5; % no. of parameters
Region(5).explain_params = {'x0','y0','z0','r_in', 'r_out'}; % explanation of parameters
Region(5).units_params={'m','m','m','m','m'};
Region(5).default=1e-6*[0 0 0 0.1 0.2];
Region(5).params=[];
Region(5).DOF=2; % allowed 2,3
%% Ellipsoid/Ellipse Shell Concentric type: ((x-x0)/a_out).^2 + ((y-y0)/b_out).^2 + ((z-z0)/c_out).^2 <= 1 AND ((x-x0)/a_in).^2 + ((y-y0)/b_in).^2 + ((z-z0)/c_in).^2 >= 1
Region(6).name='EllipsoidShellConcentric'; % shape
Region(6).no_of_params=9; % no. of parameters
Region(6).explain_params = {'x0','y0','z0','a_in','b_in','c_in','a_out','b_out','c_out'}; % explanation of parameters
Region(6).units_params={'m','m','m','m','m','m','m','m','m'};
Region(6).default=1e-6*[0 0 0 0.1 0.1 0.1 0.2 0.2 0.2];
Region(6).params=[];
Region(6).DOF=2; % allowed 2,3
%% Cube/Square Shell Concentric type a_in/2 <= abs(x-x0), abs(y-y0), abs(z-z0) <= a_out/2
Region(7).name='CubeShellConcentric'; % shape
Region(7).no_of_params=5; % no. of parameters
Region(7).explain_params = {'x0','y0','z0','a_in','a_out'}; % explanation of parameters
Region(7).units_params={'m','m','m','m','m'};
Region(7).default=1e-6*[0 0 0 0.1 0.2];
Region(7).params=[];
Region(7).DOF=2; % allowed 2,3
%% Cuboid/Rectangle Shell Concentric type a_in/2 <= abs(x-x0) <= a_out/2, b_in/2 <= abs(y-y0) <= b_out/2, c_in/2<= abs(z-z0) <= c_out/2
Region(8).name='CuboidShellConcentric'; % shape
Region(8).no_of_params=9; % no. of parameters
Region(8).explain_params = {'x0','y0','z0','a_in','b_in','c_in','a_out','b_out','c_out'}; % explanation of parameters
Region(8).units_params={'m','m','m','m','m','m','m','m','m'};
Region(8).default=1e-6*[0 0 0 0.1 0.1 0.1 0.2 0.2 0.2];
Region(8).params=[];
Region(8).DOF=2; % allowed 2,3
%% Sphere/Circle Shell Nonconcentric type: r_in.^2 <= (x-x0_in).^2 + (y-y0_in).^2 + (z-z0_in).^2 AND  (x-x0_out).^2 + (y-y0_out).^2 + (z-z0_out).^2 <= r_out.^2;
Region(9).name='SphereShellNonconcentric'; % shape
Region(9).no_of_params=8; % no. of parameters
Region(9).explain_params = {'x0_in','y0_in','z0_in','r_in', 'x0_out','y0_out','z0_out','r_out'}; % explanation of parameters
Region(9).units_params={'m','m','m','m','m','m','m','m'};
Region(9).default=1e-6*[0 0 0 0.1 0 0 0 0.2];
Region(9).params=[];
Region(9).DOF=2; % allowed 2,3
%% Ellipsoid/Ellipse Shell Nonconcentric type: ((x-x0_out)/a_out).^2 + ((y-y0_out)/b_out).^2 + ((z-z0_out)/c_out).^2 <= 1 AND ((x-x0_in)/a_in).^2 + ((y-y0_in)/b_in).^2 + ((z-z0_in)/c_in).^2 >= 1
Region(10).name='EllipsoidShellNonconcentric'; % shape
Region(10).no_of_params=12; % no. of parameters
Region(10).explain_params = {'x0_in','y0_in','z0_in','a_in','b_in','c_in','x0_out','y0_out','z0_out','a_out','b_out','c_out'}; % explanation of parameters
Region(10).units_params={'m','m','m','m','m','m','m','m','m','m','m','m'};
Region(10).default=1e-6*[0 0 0 0.1 0.1 0.1 0 0 0 0.2 0.2 0.2];
Region(10).params=[];
Region(10).DOF=2; % allowed 2,3
%% Cube/Square Shell Nonconcentric type a_in/2 <= abs(x-x0_in), abs(y-y0_in), abs(z-z0_in) AND abs(x-x0_out), abs(y-y0_out), abs(z-z0_out) <= a_out/2
Region(11).name='CubeShellNonconcentric'; % shape
Region(11).no_of_params=8; % no. of parameters
Region(11).explain_params = {'x0_in','y0_in','z0_in','a_in','x0_in','y0_in','z0_in','a_out'}; % explanation of parameters
Region(11).units_params={'m','m','m','m','m','m','m','m'};
Region(11).default=1e-6*[0 0 0 0.1 0 0 0 0.2];
Region(11).params=[];
Region(11).DOF=2; % allowed 2,3
%% Cuboid/Rectangle Shell Nonconcentric type a_in/2 <= abs(x-x0_in) AND abs(x-x0_out) <= a_out/2 AND b_in/2 <= abs(y-y0_in) AND abs(y-y0_out) <= b_out/2 AND c_in/2<= abs(z-z0_in) AND abs(z-z0_out) <= c_out/2
Region(12).name='CuboidShellNonconcentric'; % shape
Region(12).no_of_params=12; % no. of parameters
Region(12).explain_params = {'x0_in','y0_in','z0_in','a_in','b_in','c_in','x0_out','y0_out','z0_out','a_out','b_out','c_out'}; % explanation of parameters
Region(12).units_params={'m','m','m','m','m','m','m','m','m','m','m','m'};
Region(12).default=1e-6*[0 0 0 0.1 0.1 0.1 0 0 0 0.2 0.2 0.2];
Region(12).params=[];
Region(12).DOF=2; % allowed 2,3
%% Line in x-y plane y-y2=(x-x2)*(y1-y2)/(x1-x2)
Region(13).name='Line'; % shape
Region(13).no_of_params=4; % no. of parameters
Region(13).explain_params = {'x1','y1','x2','y2'}; % explanation of parameters
Region(13).units_params={'m','m','m','m'};
Region(13).default=1e-6*[-0.1 0 0.1 0];
Region(13).params=[];
Region(13).DOF=1; % allowed 1
%%

%%ARIF

%% Line in x-y plane y-y2=(x-x2)*(y1-y2)/(x1-x2)
Region(14).name='Lysosome'; % shape
Region(14).no_of_params=1; % no. of parameters
Region(14).explain_params = {'r1'}; % explanation of parameters
Region(14).units_params={'m'};
Region(14).default=1e-9*[150];
Region(14).params=[];
Region(14).DOF=1; % allowed 1

%%
%% Distribution pattern
%{
Pattern structure has the following fields:
name: 
no_of_params:
explain_params:
unit_params:normalized by min([x_max-x_min,y_max-y_min,z_max-z_min])
default:
min_dist:
%}
%{
For random patterns, we will first find many more random points than needed within the bouding box.
Then, we select the ones that are within the shape
Lastly, we select the ones that satisfy distance criterion
%}
%{
For Regular patterns, say N is the specified number of emitters
If fit_N is true,
various values of min_dist and offset_dist are used till number of emitters
is as close to N as possible
If fit_N is false,
offset_dist is set to zero and specified min_dist is used.
%}
%% Random with uniform pdf
Pattern(1).name='RandomUniform';
Pattern(1).no_of_params=0;
Pattern(1).explain_params=cell(0);
Pattern(1).units_params=cell(0);
Pattern(1).default=[];
Pattern(1).params=[];
Pattern(1).min_dist=1e-8;
%% Random with normal pdf
Pattern(2).name='RandomNormal';
Pattern(2).no_of_params=2;
Pattern(2).explain_params={'mu','sigma'};
Pattern(2).units_params={'',''};
Pattern(2).default=[0 1];
Pattern(2).params=[];
Pattern(2).min_dist=1e-8;
%% Regular with triangular gird - Bravais Lattice : http://en.wikipedia.org/wiki/Bravais_lattice
Pattern(3).name='RegularGridTriangular';
Pattern(3).no_of_params=1;
Pattern(3).explain_params={'fit_N'};
Pattern(3).units_params={'Boolean'};
Pattern(3).default=false;
Pattern(3).params=[];
Pattern(3).min_dist=1e-8;
%% Regular with square gird - Bravais Lattice : http://en.wikipedia.org/wiki/Bravais_lattice
Pattern(4).name='RegularGridSquare';
Pattern(4).no_of_params=1;
Pattern(4).explain_params={'fit_N'};
Pattern(4).units_params={'Boolean'};
Pattern(4).default=false;
Pattern(4).params=[];
Pattern(4).min_dist=1e-8;
%%
%%
%% Dynamics Library
%% Stationary
Dynamics(1).name='Stationary';
Dynamics(1).no_of_params=0;
Dynamics(1).explain_params=cell(0);
Dynamics(1).params_units=cell(0);
Dynamics(1).default=[];
Dynamics(1).params=[];
%% Flow
Dynamics(2).name='Flow';
Dynamics(2).no_of_params=6;
Dynamics(2).explain_params={'v_x','v_y','v_z','a_x','a_y','a_z'};
Dynamics(2).params_units={'m/s','m/s','m/s','m/s^2','m/s^2','m/s^2'};
Dynamics(2).default=[1e-6 1e-6 0 0 0 0 ];
Dynamics(2).params=[];
%% Diffusion
Dynamics(3).name='Diffuse';
Dynamics(3).no_of_params=1;
Dynamics(3).explain_params={'DiffuseCoef'};
Dynamics(3).params_units={'m^2/s'};
Dynamics(3).default=1e-14;
Dynamics(3).params=[];
%% DiffuseFlow
Dynamics(4).name='DiffuseFlow';
Dynamics(4).no_of_params=7;
Dynamics(4).explain_params={'DiffuseCoef','v_x','v_y','v_z','a_x','a_y','a_z'};
Dynamics(4).params_units={'m^2/s','m/s','m/s','m/s','m/s^2','m/s^2','m/s^2'};
Dynamics(4).default=[1e-14 1e-6 1e-6 0 0 0 0 ];
Dynamics(4).params=[];
%%
%% RandomWalk
Dynamics(5).name='RandomWalk';
Dynamics(5).no_of_params=3;
Dynamics(5).explain_params={'max_x','max_y','max_z'};
Dynamics(5).params_units={'nm','nm','nm'};
Dynamics(5).default=[5 5 0];
Dynamics(5).params=[];
%%

%% Random Flow
Dynamics(6).name='RandomFlow';
Dynamics(6).no_of_params=3;
Dynamics(6).explain_params={'max_x','max_y','max_z'};
Dynamics(6).params_units={'nm','nm','nm'};
Dynamics(6).default=[5 5 0];
Dynamics(6).params=[];
%%
%%
%% Flow and Stationary Mixed

%% 
Dynamics(7).name='FlowStationary';
Dynamics(7).no_of_params=6;
Dynamics(7).explain_params={'v_x','v_y','v_z','a_x','a_y','a_z'};
Dynamics(7).params_units={'m/s','m/s','m/s','m/s^2','m/s^2','m/s^2'};
Dynamics(7).default=[1e-6 1e-6 0 0 0 0 ];
Dynamics(7).params=[];
%%
%%
%% 
Dynamics(8).name='FlowRandomWalk';
Dynamics(8).no_of_params=9;
Dynamics(8).explain_params={'v_x','v_y','v_z','a_x','a_y','a_z','max_x','max_y','max_z'};
Dynamics(8).params_units={'m/s','m/s','m/s','m/s^2','m/s^2','m/s^2','nm','nm','nm'};
Dynamics(8).default=[1e-6 1e-6 0 0 0 0 5 5 0];
Dynamics(8).params=[];
%%
%%
%% 
Dynamics(9).name='RandomWalkCurve';
Dynamics(9).no_of_params=5;
Dynamics(9).explain_params={'max_x','max_y','max_z','a','b'};
Dynamics(9).params_units={'nm','nm','nm','c','c'};
Dynamics(9).default=[5 5 0 2 3];
Dynamics(9).params=[];
%%
%%
%% RandomWalk in a Circle
Dynamics(10).name='RandomWalkCircle';
Dynamics(10).no_of_params=1;
Dynamics(10).explain_params={'max_r'};
Dynamics(10).params_units={'nm'};
Dynamics(10).default=[5];
Dynamics(10).params=[];
%%
%%
%%
%% CircularMotion
Dynamics(11).name='CircularMotion';
Dynamics(11).no_of_params=1;
Dynamics(11).explain_params={'max_r'};
Dynamics(11).params_units={'nm'};
Dynamics(11).default=[5];
Dynamics(11).params=[];
%%

%%
%% CircularMotion
Dynamics(12).name='CircularMotionRange';
Dynamics(12).no_of_params=2;
Dynamics(12).explain_params={'min_r','max_r'};
Dynamics(12).params_units={'nm','nm'};
Dynamics(12).default=[5 10];
Dynamics(12).params=[];
%%

%% RandomWalk inside a Circle
Dynamics(13).name='RandomWalkInsideCircle';
Dynamics(13).no_of_params=2;
Dynamics(13).explain_params={'max_v','max_r'};
Dynamics(13).params_units={'nm','nm'};
Dynamics(13).default=[5 15];
Dynamics(13).params=[];
%%

%% RandomWalk Small Region
Dynamics(14).name='RandomWalkSmallregion';
Dynamics(14).no_of_params=2;
Dynamics(14).explain_params={'max_v','max_r'};
Dynamics(14).params_units={'nm','nm'};
Dynamics(14).default=[5 5];
Dynamics(14).params=[];
%%
%%
%% Blinking Library
%% Binomial
Blinking(1).name='Binomial';
Blinking(1).no_of_params=4;
Blinking(1).explain_params={'photonspersec', 'p_on', 'tau_blink', 'p_bleach'};
Blinking(1).units_params={'/second' '','seconds',''};
Blinking(1).default=[1e5 0.3 1e-3 0.1];
Blinking(1).params=[];
%% Hidden Markov
Blinking(2).name='HiddenMarkov';
Blinking(2).no_of_params=6;
Blinking(2).explain_params={'photonspersec', 'p_on','p_off','p_on_off','p_off_on','p_bleach'};
Blinking(2).units_params={'/second' '','','','',''};
Blinking(2).default=[1e5 0.3 0.3 0.3 0.3 0.3];
Blinking(2).params=[];
%% Poisson - to modify later
Blinking(3).name='Poisson';
Blinking(3).no_of_params=6;
Blinking(3).explain_params={'photonspersec', 'p_on' , 'tau_blink', 'p_bleach'};
Blinking(3).units_params={'/second' '','seconds',''};
Blinking(3).default=[1e5 0.3 1e-3 0.1];
Blinking(3).params=[];
%% Poisson - to modify later
Blinking(4).name='Poisson_ExpBleach';
Blinking(4).no_of_params=6;
Blinking(4).explain_params={'photonspersec', 'p_on' , 'tau_blink', 'tau_bleach'};
Blinking(4).units_params={'/second' '','seconds','seconds'};
Blinking(4).default=[1e5 0.3 1e-3 0.01];
Blinking(4).params=[];
%%
%%
%% Save libraries
save('Libraries','Region','Pattern','Dynamics','Blinking');