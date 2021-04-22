function [coords,emission,NoOfEmitters,ExceptionString]=...
    generate_sample_and_dynamics(Region,Pattern,Dynamics,Blinking,NoOfEmitters,time_step,ObservationTime,xyz_span)
%% Generate initial state
if isempty(Region.params)
    Region.params=Region.default;
end
for param_no=1:length(Region.explain_params)
    %%% get parameters
    eval([Region.explain_params{param_no} '=' num2str(Region.params(param_no)) ';']);
end
if isempty(Pattern.params)
    Pattern.params=Pattern.default;
end
for param_no=1:length(Pattern.explain_params)
    %%% get parameters
    eval([Pattern.explain_params{param_no} '=' num2str(Pattern.params(param_no)) ';']);
end
ExceptionString='';
switch upper(Region.name)
    %%% Sphere
    case upper('Sphere')
        switch Region.DOF
            % Sphere edge
            case 1
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                        x=[];y=[];z=[];force_N=true;
                        while(force_N)
                            th=random('Uniform',0,2*pi,[10*NoOfEmitters,1]);
                            x_t = r.*cos(th) + x0;
                            y_t = r.*sin(th) + y0;
                            z_t = zeros(size(th));
                            x=[x(:); x_t(:)];y=[y(:); y_t(:)];z=[z(:); z_t(:)];
                            min_dist=min(generate_dist([x y z],[x y z]));
                            keep=find(min_dist)>Pattern.min_dist;
                            if length(keep)>=NoOfEmitters
                                force_N=false;
                                x=x(keep(1:NoOfEmitters));y=y(keep(1:NoOfEmitters));z=z(keep(1:NoOfEmitters));
                            end
                        end
                    case upper('RandomNormal')
                        disp([Region.name ', DOF = ' num2str(Region.DOF) ', ' Pattern.name]);
                        disp('Generating as RandomUniform');
                        ExceptionString=[ExceptionString 'Generating as RandomUniform;'];
                        x=[];y=[];z=[];force_N=true;
                        while(force_N)
                            th=random('Uniform',0,2*pi,[10*NoOfEmitters,1]);
                            x_t = r.*cos(th) + x0;
                            y_t = r.*sin(th) + y0;
                            z_t = zeros(size(th));
                            x=[x(:); x_t(:)];y=[y(:); y_t(:)];z=[z(:); z_t(:)];
                            min_dist=min(generate_dist([x y z],[x y z]));
                            keep=find(min_dist)>Pattern.min_dist;
                            if length(keep)>=NoOfEmitters
                                force_N=false;
                                x=x(keep(1:NoOfEmitters));y=y(keep(1:NoOfEmitters));z=z(keep(1:NoOfEmitters));
                            end
                        end
                    case upper('RegularGridTriangular')
                        disp([Region.name ', DOF = ' num2str(Region.DOF) ', ' Pattern.name]);
                        disp('Generating as RegularAngular');
                        ExceptionString=[ExceptionString 'Generating as RegularAngular;'];
                        if (isempty(Pattern.params) && Pattern.default==true) || Pattern.params==true
                            th=linspace(2*pi/NoOfEmitters,2*pi,NoOfEmitters);
                        else
                            d_th=r/Pattern.min_dist;
                            NoOfEmitters=floor(2*pi/d_th);
                            th=linspace(2*pi/NoOfEmitters,2*pi,NoOfEmitters);
                            ExceptionString=[ExceptionString 'NoOfEmitters updated for min_dist'];
                        end
                        x = r*cos(th) + x0;
                        y = r*sin(th) + y0;
                        z = zeros(size(th));
                    case upper('RegularGridSquare')
                        disp([Region.name ', DOF = ' num2str(Region.DOF) ', ' Pattern.name]);
                        disp('Generating as RegularAngular');
                        ExceptionString=[ExceptionString 'Generating as RegularAngular;'];
                        if (isempty(Pattern.params) && Pattern.default==true) || Pattern.params==true
                            th=linspace(2*pi/NoOfEmitters,2*pi,NoOfEmitters);
                        else
                            d_th=r/Pattern.min_dist;
                            NoOfEmitters=floor(2*pi/d_th);
                            th=linspace(2*pi/NoOfEmitters,2*pi,NoOfEmitters);
                            ExceptionString=[ExceptionString 'NoOfEmitters updated for min_dist'];
                        end
                        x = r*cos(th) + x0;
                        y = r*sin(th) + y0;
                        z = zeros(size(th));
                end
                
                % Sphere surface
            case 2
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                        x=[];y=[];z=[];force_N=true;
                        while(force_N)
                            th=random('Uniform',0,2*pi,[10*NoOfEmitters,1]);
                            r_em=random('Uniform',0,r,[10*NoOfEmitters,1]);
                            x_t = r_em.*cos(th) + x0;
                            y_t = r_em.*sin(th) + y0;
                            z_t = zeros(size(th));
                            x=[x(:); x_t(:)];y=[y(:); y_t(:)];z=[z(:); z_t(:)];
                            min_dist=min(generate_dist([x y z],[x y z]));
                            keep=find(min_dist)>Pattern.min_dist;
                            if length(keep)>=NoOfEmitters
                                force_N=false;
                                x=x(keep(1:NoOfEmitters));y=y(keep(1:NoOfEmitters));z=z(keep(1:NoOfEmitters));
                            end
                        end
                    case upper('RandomNormal')
                        x=[];y=[];z=[];force_N=true;
                        while(force_N)
                            th=random('Uniform',0,2*pi,[10*NoOfEmitters,1]);
                            r_em=abs(random('Normal',Pattern.params(1)*r,Pattern.params(2)*r,[10*NoOfEmitters,1]));
                            x_t = r_em.*cos(th) + x0;
                            y_t = r_em.*sin(th) + y0;
                            z_t = zeros(size(th));
                            x=[x(:); x_t(:)];y=[y(:); y_t(:)];z=[z(:); z_t(:)];
                            min_dist=min(generate_dist([x y z],[x y z]));
                            keep=find(min_dist)>Pattern.min_dist;
                            if length(keep)>=NoOfEmitters
                                force_N=false;
                                x=x(keep(1:NoOfEmitters));y=y(keep(1:NoOfEmitters));z=z(keep(1:NoOfEmitters));
                            end
                        end
                    case upper('RegularGridTriangular')
                        force_N=(isempty(Pattern.params) && Pattern.default==true) || Pattern.params==true;
                        [x,y,z,ExceptionString] = function_generate_regular_grid(force_N,NoOfEmitters,Region.DOF,2,2,Pattern.min_dist,'ellipsoid_shell',[x0 y0],r,[],[],ExceptionString);
                    case upper('RegularGridSquare')
                        force_N=(isempty(Pattern.params) && Pattern.default==true) || Pattern.params==true;
                        [x,y,z,ExceptionString] = function_generate_regular_grid(force_N,NoOfEmitters,Region.DOF,2,1,Pattern.min_dist,'ellipsoid_shell',[x0 y0],r,[],[],ExceptionString);
                end
                % Sphere volume
            case 3
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                end
            otherwise
        end
    case upper('Ellipsoid')
        switch Region.DOF
            case 1
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 2
                switch upper(Pattern.name)
                    case ('RandomUniform')
                    case ('RandomNormal')
                    case ('RegularGridTriangular')
                    case ('RegularGridSquare')
                    otherwise
                end
            case 3
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            otherwise
        end
    case upper('Cube')
        switch Region.DOF
            case 1
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 2
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 3
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            otherwise
        end
    case upper('Cuboid')
        switch Region.DOF
            case 1
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 2
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 3
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            otherwise
        end
    case upper('SphereShellConcentric')
        switch Region.DOF
            case 1
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 2
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 3
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            otherwise
        end
    case upper('EllipsoidShellConcentric')
        switch Region.DOF
            case 1
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 2
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 3
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            otherwise
        end
    case upper('CubeShellConcentric')
        switch Region.DOF
            case 1
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 2
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 3
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            otherwise
        end
    case upper('CuboidShellConcentric')
        switch Region.DOF
            case 1
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 2
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 3
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            otherwise
        end
    case upper('SphereShellNonconcentric')
        switch Region.DOF
            case 1
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 2
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 3
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            otherwise
        end
    case upper('EllipsoidShellNonconcentric')
        switch Region.DOF
            case 1
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 2
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 3
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            otherwise
        end
    case upper('CubeShellNonconcentric')
        switch Region.DOF
            case 1
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 2
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 3
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            otherwise
        end
    case upper('CuboidShellNonconcentric')
        switch Region.DOF
            case 1
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 2
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 3
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            otherwise
        end
    case upper('Line')
        switch Region.DOF
            case 1
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 2
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            case 3
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                    case upper('RandomNormal')
                    case upper('RegularGridTriangular')
                    case upper('RegularGridSquare')
                    otherwise
                end
            otherwise
        end
    otherwise
end
coords=[x(:) y(:) z(:)];
NoOfEmitters=size(coords,1);
%% Emitter's kinetics dynamics
if isempty(Dynamics.params)
    Dynamics.params=Dynamics.default;
end
for param_no=1:length(Dynamics.explain_params)
    %%% get parameters
    eval([Dynamics.explain_params{param_no} '=' num2str(Dynamics.params(param_no)) ';']);
end
switch upper(Dynamics.name)
    case upper('stationary')
        for t=time_step:time_step:ObservationTime
            coords(:,:,end+1)=coords(:,:,end);
        end
    case upper('Flow')
        for t=time_step:time_step:ObservationTime
            coords(:,:,end+1)=coords(:,:,end) + repmat([v_x v_y v_z]*time_step,size(coords,1),1);
        end
    case upper('Diffuse')
        for t=time_step:time_step:ObservationTime
            perturb=random('Normal',0,sqrt(2*DiffuseCoef*t),[size(coords,1) size(coords,2)]);
            if isempty(find(coords(:,3,end)~=0))
                perturb(:,3)=0;
            end
            particles_temp=coords(:,:,end)+perturb;
            need_check=true;
            while(need_check)
                indices_to_change=find((abs(particles_temp(:,1))>xyz_span(1)/2) | (abs(particles_temp(:,2))>xyz_span(2)/2) | (abs(particles_temp(:,3))>xyz_span(3)/2));
                if isempty(indices_to_change)
                    need_check=false;
                else
                    perturb(indices_to_change,:)=random('Normal',0,sqrt(2*DiffuseCoef*t),[length(indices_to_change) size(coords,2)]);
                    if isempty(find(coords(:,3,end)~=0))
                        perturb(:,3)=0;
                    end
                    particles_temp=coords(:,:,end)+perturb;
                end
            end
            coords(:,:,end+1)=coords(:,:,end) + perturb;
        end
    case upper('DiffuseFlow')
        for t=time_step:time_step:ObservationTime
            perturb=random('Normal',0,sqrt(2*DiffuseCoef*t),[size(coords,1) size(coords,2) size(coords,3)]);
            if isempty(find(coords(:,3,end)~=0))
                perturb(:,3)=0;
            end
            particles_temp=coords(:,:,end)+perturb;
            need_check=true;
            while(need_check)
                indices_to_change=find((abs(particles_temp(:,1))>xyz_span(1)/2) | (abs(particles_temp(:,2))>xyz_span(2)/2) | (abs(particles_temp(:,3))>xyz_span(3)/2));
                if isempty(indices_to_change)
                    need_check=false;
                else
                    perturb=random('Normal',0,sqrt(2*DiffuseCoef*t),[size(coords,1) size(coords,2) size(coords,3)]);
                    particles_temp=coords(:,:,end)+perturb;
                end
            end
            perturbflow=repmat([v_x v_y v_z]*time_step,size(coords,1),1);
            coords(:,:,end+1)=coords(:,:,end) + perturb+perturbflow;
        end
end
%% blinking state
if isempty(Blinking.params)
    Blinking.params=Blinking.default;
end
for param_no=1:length(Blinking.explain_params)
    %%% get parameters
    eval([Blinking.explain_params{param_no} '=' num2str(Blinking.params(param_no)) ';']);
end
NoOfEmitters=size(coords,1);
NoOfTimeSteps=size(coords,3);
photonsperstep=photonspersec*time_step;
T_chunk=round(tau_blink/time_step);

emission=abs(round(random('Normal',photonsperstep,sqrt(photonsperstep),[NoOfEmitters NoOfTimeSteps])));
moment_bleach=NoOfTimeSteps-random('geo',p_bleach,[NoOfEmitters,1]);
no_of_ONs=random('Binomial',T_chunk,p_on,[NoOfEmitters floor(NoOfTimeSteps/T_chunk)]);
blinking=zeros([NoOfEmitters NoOfTimeSteps]);
for n=1:NoOfEmitters
    for nn=1:floor(NoOfTimeSteps/T_chunk)
        if no_of_ONs(n,nn)~=0
            if no_of_ONs(n,nn)~=T_chunk
                tt=randperm(T_chunk-no_of_ONs(n,nn),1);
            else
                tt=0;
            end
            blinking(n,(tt+(nn-1)*T_chunk+1):(tt+(nn-1)*T_chunk+no_of_ONs(nn)))=1;
        end
    end
    blinking(n,moment_bleach(n):end)=0;
end
emission=emission.*blinking(:,1:NoOfTimeSteps);
end

function [x,y,z,ExceptionString] = function_generate_regular_grid(force_N,NoOfEmitters,DOF,dim,type,min_dist,shape,centers_out,shape_dim_out,centers_in,shape_dim_in,ExceptionString)
if strcmpi(shape,'ellipsoid_shell')
    span=2*shape_dim_out;
else
    span=shape_dim_out;
end
if length(span)==1
    span=span*[1 1 1];
end
if length(centers_out)<3
    centers_out(3)=0;
end
[x,y,z] = function_Bravais_grid(dim,type,min_dist,span);
x=x+centers_out(1);y=y+centers_out(2);z=z+centers_out(3);
keep=eval(['function_check_' shape '([x(:),y(:),z(:)],centers_out,shape_dim_out,centers_in,shape_dim_in);']);
if force_N
    ExceptionString=[ExceptionString 'min_dist Changed;'];
    N_current=length(keep);
    if N_current~=NoOfEmitters
        fit_N=true;
        d_current=min_dist;
    else
        fit_N=false;
    end
    while (fit_N)
        d_prev=d_current;
        N_previous=N_current;
        d_current=min_dist*(N_current/NoOfEmitters)^(1/DOF);
        [x,y,z] = function_Bravais_grid(dim,type,d_current,span);
        x=x+centers_out(1);y=y+centers_out(2);z=z+centers_out(3);
        keep=eval(['function_check_' shape '([x(:),y(:),z(:)],centers_out,shape_dim_out,centers_in,shape_dim_in);']);
        N_current=length(keep);
        fit_N=false;
        if N_previous==N_current
            fit_N=false;
        end
    end
else
    ExceptionString=[ExceptionString 'NoOfEmitters Changed;'];
end
x=x(keep);
y=y(keep);
z=z(keep);
end
function [x,y,z] = function_Bravais_grid(dim,type,dist,span)
% dim: 3D or 3D
% type: 1: cartesian (square or cubic), 2: hexagonal
% dist: distance of unit cell
% span: [span_x,span_y,span_z]
switch dim
    case 2
        switch type
            case 1
                a_1=dist*[1 0 0];
                a_2=dist*[0 1 0];
            case 2
                a_1=dist*[1 0 0];
                a_2=dist*[cosd(120) sind(120) 0];
        end
        n_x=ceil((span(1)/2)/dist);
        n_y=ceil((span(2)/2)/dist);
        [nn1,nn2]=meshgrid(-n_x:n_x,-n_y:n_y);
        xy=nn1(:)*a_1 + nn2(:) *a_2;
        x=xy(:,1);
        y=xy(:,2);
        z=xy(:,3);
    case 3
        switch type
            case 1
                a_1=dist*[1 0 0];
                a_2=dist*[0 1 0];
                a_3=dist*[0 0 1];
            case 2
                a_1=dist*[1 0 0];
                a_2=dist*[cosd(120) sind(120) 0];
                a_3=dist*[0 0 1];
        end
        n_x=ceil((span(1)/2)/dist);
        n_y=ceil((span(2)/2)/dist);
        n_z=ceil((span(2)/2)/dist);
        [nn1,nn2,nn3]=meshgrid(-n_x:n_x,-n_y:n_y,-n_z:n_z);
        xy=nn1(:)*a_1 + nn2(:)*a_2 + nn3(:)*a_3;
        x=xy(:,1);
        y=xy(:,2);
        z=xy(:,3);
end
end
function inside=function_check_ellipsoid_shell(coords,centers_out,semiaxes_out,centers_in,semiaxes_in)
x=coords(:,1);
y=coords(:,2);
if size(coords,2)==3
    z=coords(:,3);
else
    z=zeros(size(x));
    c_out=0;c_in=0;z0_out=0;z0_in=0;
end
% centers_out
if isempty(centers_out) || length(centers_out)==1
    x0_out=0;y0_out=0;z0_out=0;
elseif length(centers_out)==2
    x0_out=centers_out(1);y0_out=centers_out(2);z0_out=0;
else
    x0_out=centers_out(1);y0_out=centers_out(2);z0_out=centers_out(3);
end
% semiaxes_out
if length(semiaxes_out)==1
    a_out=semiaxes_out;b_out=semiaxes_out;c_out=semiaxes_out;
elseif length(semiaxes_out)==2
    a_out=semiaxes_out(1);b_out=semiaxes_out(2);c_out=0;
else
    a_out=semiaxes_out(1);b_out=semiaxes_out(2);c_out=semiaxes_out(3);
end
% centers_in
if isempty(centers_in) || length(centers_in)==1
    x0_in=0;y0_in=0;z0_in=0;
elseif length(centers_in)==2
    x0_in=centers_in(1);y0_in=centers_in(2);z0_in=0;
else
    x0_in=centers_in(1);y0_in=centers_in(2);z0_in=centers_in(3);
end
% semiaxes_in
if isempty(semiaxes_in)
    a_in=0;b_in=0;c_in=0;
elseif length(semiaxes_in)==1
    a_in=semiaxes_in;b_in=semiaxes_in;c_in=semiaxes_in;
elseif length(semiaxes_in)==2
    a_in=semiaxes_in(1);b_in=semiaxes_in(2);c_in=0;
else
    a_in=semiaxes_in(1);b_in=semiaxes_in(2);c_in=semiaxes_in(3);
end
if size(coords(2))==2
    in_check=((x-x0_in)/a_in).^2 + ((y-y0_in)/b_in).^2;
    in_check(isnan(in_check))=1;
else
    in_check=((x-x0_in)/a_in).^2 + ((y-y0_in)/b_in).^2 + ((z-z0_in)/c_in).^2;
    in_check(isnan(in_check))=1;
end
if size(coords(2))==2
    out_check=((x-x0_out)/a_out).^2 + ((y-y0_out)/b_out).^2;
else
    out_check=((x-x0_out)/a_out).^2 + ((y-y0_out)/b_out).^2 + ((z-z0_out)/c_out).^2;
end
inside=find(in_check>=1 & out_check<=1);
end
function inside=function_check_cuboid_shell(coords,centers_out,sides_out,centers_in,sides_in)
x=coords(:,1);
y=coords(:,2);
if size(coords,2)==3
    z=coords(:,3);
else
    z=zeros(size(x));
    c_out=0;c_in=0;z0_out=0;z0_in=0;
end
% centers_out
if isempty(centers_out) || length(centers_out)==1
    x0_out=0;y0_out=0;z0_out=0;
elseif length(centers_out)==2
    x0_out=centers_out(1);y0_out=centers_out(2);z0_out=0;
else
    x0_out=centers_out(1);y0_out=centers_out(2);z0_out=centers_out(3);
end
% sides_out
if length(sides_out)==1
    a_out=sides_out;b_out=sides_out;c_out=sides_out;
elseif length(sides_out)==2
    a_out=sides_out(1);b_out=sides_out(2);c_out=0;
else
    a_out=sides_out(1);b_out=sides_out(2);c_out=sides_out(3);
end
% centers_in
if isempty(centers_in) || length(centers_in)==1
    x0_in=0;y0_in=0;z0_in=0;
elseif length(centers_in)==2
    x0_in=centers_in(1);y0_in=centers_in(2);z0_in=0;
else
    x0_in=centers_in(1);y0_in=centers_in(2);z0_in=centers_in(3);
end
% sides_in
if isempty(sides_in)
    a_in=0;b_in=0;c_in=0;
elseif length(sides_in)==1
    a_in=sides_in;b_in=sides_in;c_in=sides_in;
elseif length(sides_in)==2
    a_in=sides_in(1);b_in=sides_in(2);c_in=0;
else
    a_in=sides_in(1);b_in=sides_in(2);c_in=sides_in(3);
end
if size(coords(2))==2
    in_check=(abs(x-x_in)>=(a_in)/2) | (abs(y-y_in)>=(b_in)/2);
else
    in_check=(abs(x-x_in)>=(a_in)/2) | (abs(y-y_in)>=(b_in)/2)  | (abs(z-z_in)>=(c_in)/2);
end
if size(coords(2))==2
    out_check=(abs(x-x_out)<=(a_out)/2) | (abs(y-y_out)<=(b_out)/2);
else
    out_check=(abs(x-x_out)<=(a_out)/2) | (abs(y-y_out)<=(b_out)/2)  | (abs(z-z_out)<=(c_out)/2);
end
inside=find(in_check & out_check);
end
function dist=generate_dist(coords1,coords2)
x1=coords1(:,1);
y1=coords1(:,2);
if size(coords1,2)==3
    z1=coords1(:,3);
else
    z1=zeros(size(x1));
end
x2=coords2(:,1);
y2=coords2(:,2);
if size(coords2,2)==3
    z2=coords2(:,3);
else
    z2=zeros(size(x2));
end
[xx1,xx2]=meshgrid(x1,x2);[yy1,yy2]=meshgrid(y1,y2);[zz1,zz2]=meshgrid(z1,z2);
dist=sqrt((xx1-xx2).^2 + (yy1-yy2).^2 + (zz1-zz2).^2);
dist(dist==0)=Inf;
end