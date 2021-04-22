function [coords,emission,NoOfEmitters,ExceptionString]=...
    generate_sample_and_dynamics(Region,Pattern,Dynamics,Blinking,...
    NoOfEmitters,time_step,ObservationTime,xyz_span,flag_emission_same)
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
                            d_th=Pattern.min_dist/r;
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
                            d_th=Pattern.min_dist/r;
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
                        x=[];y=[];z=[];force_N=true;
                        while(force_N)
                            th=random('Uniform',0,2*pi,[10*NoOfEmitters,1]);
                            x_t = a.*cos(th) + x0;
                            y_t = b.*sin(th) + y0;
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
                            x_t = a.*cos(th) + x0;
                            y_t = b.*sin(th) + y0;
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
                            d_th=b/Pattern.min_dist;
                            NoOfEmitters=floor(2*pi/d_th);
                            th=linspace(2*pi/NoOfEmitters,2*pi,NoOfEmitters);
                            ExceptionString=[ExceptionString 'NoOfEmitters updated for min_dist'];
                        end
                        x = a*cos(th) + x0;
                        y = b*sin(th) + y0;
                        z = zeros(size(th));
                    case upper('RegularGridSquare')
                        disp([Region.name ', DOF = ' num2str(Region.DOF) ', ' Pattern.name]);
                        disp('Generating as RegularAngular');
                        ExceptionString=[ExceptionString 'Generating as RegularAngular;'];
                        if (isempty(Pattern.params) && Pattern.default==true) || Pattern.params==true
                            th=linspace(2*pi/NoOfEmitters,2*pi,NoOfEmitters);
                        else
                            d_th=b/Pattern.min_dist;
                            NoOfEmitters=floor(2*pi/d_th);
                            th=linspace(2*pi/NoOfEmitters,2*pi,NoOfEmitters);
                            ExceptionString=[ExceptionString 'NoOfEmitters updated for min_dist'];
                        end
                        x = a*cos(th) + x0;
                        y = b*sin(th) + y0;
                        z = zeros(size(th));
                    otherwise
                end
            case 2
                switch upper(Pattern.name)
                    case ('RandomUniform')
                        x=[];y=[];z=[];force_N=true;
                        while(force_N)
                            th=random('Uniform',0,2*pi,[10*NoOfEmitters,1]);
                            a_em=random('Uniform',0,a,[10*NoOfEmitters,1]);
                            b_em=random('Uniform',0,b,[10*NoOfEmitters,1]);
                            x_t = a_em.*cos(th) + x0;
                            y_t = b_em.*sin(th) + y0;
                            z_t = zeros(size(th));
                            x=[x(:); x_t(:)];y=[y(:); y_t(:)];z=[z(:); z_t(:)];
                            min_dist=min(generate_dist([x y z],[x y z]));
                            keep=find(min_dist)>Pattern.min_dist;
                            if length(keep)>=NoOfEmitters
                                force_N=false;
                                x=x(keep(1:NoOfEmitters));y=y(keep(1:NoOfEmitters));z=z(keep(1:NoOfEmitters));
                            end
                        end
                    case ('RandomNormal')
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
                x=linspace(x1,x2,NoOfEmitters);
                y=linspace(y1,y2,NoOfEmitters);
                z=zeros(size(x));
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                        theta_slope=atan2((y2-y1),(x2-x1));
                        r=sqrt((y2-y1)^2+(x2-x1)^2)*rand(1,NoOfEmitters);
                        x=r*cos(theta_slope)+x1;
                        y=r*sin(theta_slope)+y1;
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
        
        
        %Arif Lysosome
        case upper('Lysosome')
        switch Region.DOF
            
            case 1
                Nmax = NoOfEmitters;
                R = r1;
                x0 = 0; % Center of the circle in the x direction.
                y0 = 0; % Center of the circle in the y direction.
                        % Now create the set of points.
                %t = 2*pi*rand(n,1);
                %r = R*sqrt(rand(n,1));
                %x = x0 + r.*cos(t);
                %y = y0 + r.*sin(t);
                %z=zeros(size(x));
                for n=1:Nmax
                    r2(n) = R*sqrt(rand(1,1));
                    theta2(n)=2*pi*rand(1,1);
                    x(n)=x0+r2(n)*cos(theta2(n));
                    y(n)=y0+r2(n)*sin(theta2(n));
                end
                z=zeros(size(x));
                %x=linspace(x1,x2,NoOfEmitters);
                %y=linspace(y1,y2,NoOfEmitters);
                %z=zeros(size(x));
                switch upper(Pattern.name)
                    case upper('RandomUniform')
                        n = NoOfEmitters;
                        R = r1;
                        x0 = 0; % Center of the circle in the x direction.
                        y0 = 0; % Center of the circle in the y direction.
                         % Now create the set of points.
                        for n=1:Nmax
                            r2(n) = R*sqrt(rand(1,1));
                            theta2(n)=2*pi*rand(1,1);
                            x(n)=x0+r2(n)*cos(theta2(n));
                            y(n)=y0+r2(n)*sin(theta2(n));
                        end
                        z=zeros(size(x));
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
%         for t=time_step:time_step:ObservationTime
%             coords(:,:,end+1)=coords(:,:,end);
%         end
        pos_x=randi([-250 250]);
        pos_y=randi([-250 250]);
        coords(:,:,end)=coords(:,:,end) + repmat([(1e-9)*pos_x (1e-9)*pos_y 0],size(coords,1),1);
        tt=repmat(coords(:,:,end),[1 1 length(time_step:time_step:ObservationTime)+1]);
        coords(:,:,end:end+(size(tt,3)-1))=tt;
    case upper('Flow')
        for t=time_step:time_step:ObservationTime
            coords(:,:,end+1)=coords(:,:,end) + repmat([v_x v_y v_z]*time_step,size(coords,1),1);
        end
    case upper('FlowStationary')
        event_time=length(time_step:time_step:ObservationTime);
        event_time=event_time/2;
        counter=1;
        for t=time_step:time_step:ObservationTime
            if(counter<event_time)
                coords(:,:,end+1)=coords(:,:,end) + repmat([v_x v_y v_z]*time_step,size(coords,1),1);
            else
                coords(:,:,end+1)=coords(:,:,end);
            end
            counter=counter+1;
        end
    case upper('FlowRandomWalk')
        event_time=length(time_step:time_step:ObservationTime);
        event_time=event_time/2;
        counter=1;
        for t=time_step:time_step:ObservationTime
            if(counter<event_time)
                coords(:,:,end+1)=coords(:,:,end) + repmat([v_x v_y v_z]*time_step,size(coords,1),1);
            else
                m_x=randi([-max_x,max_x]);
                m_y=randi([-max_y,max_y]);
                if(max_z>0)
                    m_z=randi([-max_z,max_z]);
                else
                    m_z=0;
                end
                coords(:,:,end+1)=coords(:,:,end) + repmat([(1e-9)*m_x (1e-9)*m_y (1e-9)*m_z],size(coords,1),1);
 
            end
            counter=counter+1;
        end
    case upper('RandomFlow')
        pos_x=randi([-250 250]);
        pos_y=randi([-250 250]);
        coords(:,:,end)=coords(:,:,end) + repmat([(1e-9)*pos_x (1e-9)*pos_y 0],size(coords,1),1);
        for t=time_step:time_step:ObservationTime
            
            m_x=randi([-max_x,max_x]);
            m_y=randi([-max_y,max_y]);
            if(max_z>0)
                m_z=randi(max_z);
            else
                m_z=0;
            end
            coords(:,:,end+1)=coords(:,:,end) + repmat([(1e-9)*m_x (1e-9)*m_y (1e-9)*m_z],size(coords,1),1);
        end
    case upper('RandomWalkCircle')
        for t=time_step:time_step:ObservationTime
            
            m_x=randi([-max_r,max_r]);
            m_y=randi([-max_r,max_r]);
            m_z=0;
            coords(:,:,end+1)=coords(:,:,end)-coords(:,:,end) + repmat([(1e-9)*m_x (1e-9)*m_y (1e-9)*m_z],size(coords,1),1);
        end
        Params1.max_r=m_x;
        Params1.min_r=m_y;
        
    case upper('RandomWalkInsideCircle')
        %Random initial position
        pos_x=randi([0 250]);
        pos_y=randi([0 250]);
        m_z=0;
        cx=mean(coords(:,1,end));
        cy=mean(coords(:,2,end));
        for t=time_step:time_step:ObservationTime
           r2 = max_r*sqrt(rand(1,1))*(1e-9);
           theta2=2*pi*rand(1,1);
           m_x=cx+r2*cos(theta2);
           m_y=cy+r2*sin(theta2);
           
           dis_x=cx-m_x;
           dis_y=cy-m_y;
           coords(:,:,end+1)=coords(:,:,1)+repmat([dis_x dis_y (1e-9)*m_z],size(coords,1),1);
           
        end
        Params1.max_r=m_x;
        Params1.min_r=m_y;
        
       case upper('RandomWalkSmallRegion')
        %Random initial position
        pos_x=randi([0 250]);
        pos_y=randi([0 250]);
        m_z=0;
        cx=mean(coords(:,1,end));
        cy=mean(coords(:,2,end));
        for t=time_step:time_step:ObservationTime
            do=true;
            while do
                cx1=mean(coords(:,1,end));
                cy1=mean(coords(:,2,end));
                %del_v = max_v*sqrt(rand(1,1));
                %del_v
                dis_x=randi([-max_v max_v])*(1e-9);
                dis_y=randi([-max_v max_v])*(1e-9);
                cx2=cx1+dis_x*(1e-9);
                cy2=cy1+dis_y*(1e-9);
                dx=cx-cx2;
                dy=cy-cy2;
                dis=dx^2+dy^2;
                if(dis<(max_r^2)*(1e-9))
                    do=false;
                end
                
            end
           coords(:,:,end+1)=coords(:,:,1)+repmat([dis_x dis_y (1e-9)*m_z],size(coords,1),1);
           
        end
%        Params1.max_r=m_x;
%        Params1.min_r=m_y;
    case upper('CircularMotion')
        a=0;
        cx=mean(coords(:,1,end));
        cy=mean(coords(:,2,end));
        %Random shift
        rx=randi([1 500]);
        ry=randi([1 500]);
        coords(:,:,end)=coords(:,:,end)+repmat([(1e-9)*rx (1e-9)*ry (1e-9)*0],size(coords,1),1);
        for t=time_step:time_step:ObservationTime
            
            cx=mean(coords(:,1,end));
            cy=mean(coords(:,2,end));
            m_x=cos(a)*max_r;
            m_y=sin(a)*max_r;
            m_z=0;
            dis_x=cx-m_x;
            dis_y=cy-m_y;
            coords(:,:,end+1)=coords(:,:,end)+repmat([(1e-9)*dis_x (1e-9)*dis_y (1e-9)*m_z],size(coords,1),1);
            a=a+0.001*randi([1 20]);
            if(a>=2*pi)
                a=0;
            end
        end
        Params1.max_r=max_r;
        
    case upper('CircularMotionRange')
        pos_x=randi([0 250]);
        pos_y=randi([0 250]);
        pos_z=0;
        a=0;
        cx=mean(coords(:,1,end));
        cy=mean(coords(:,2,end));
        
        for t=time_step:time_step:ObservationTime
            cx=mean(coords(:,1,end));
            cy=mean(coords(:,2,end));
            r=randi([min_r max_r])/100;
            m_x=cos(a)*r;
            m_y=sin(a)*r;
            m_z=0;
            dis_x=cx-m_x;
            dis_y=cy-m_y;
            %coords(:,:,end+1)=coords(:,:,end)-coords(:,:,end) +repmat([(1e-9)*pos_x (1e-9)*pos_y (1e-9)*pos_z],size(coords,1),1)+ repmat([(1e-9)*m_x (1e-9)*m_y (1e-9)*m_z],size(coords,1),1);
            coords(:,:,end+1)=coords(:,:,end)+repmat([(1e-9)*dis_x (1e-9)*dis_y (1e-9)*m_z],size(coords,1),1);
            
             a=a+0.001*randi([1 20]);
            if(a>=2*pi)
                a=0;
            end
        end
        Params1.max_r=max_r;
        Params1.min_r=min_r;
    case upper('LinearRandomWalk')
        for t=time_step:time_step:ObservationTime
            
            m_x=randi(max_x);
            m_y=randi(max_y);
            if(max_z>0)
                m_z=randi(max_z);
            else
                m_z=0;
            end
            coords(:,:,end+1)=coords(:,:,end) + repmat([(1e-9)*m_x (1e-9)*m_y (1e-9)*m_z],size(coords,1),1);
        end
    case upper('RandomWalkCurve')
        pos_x=randi([-250 250]);
        pos_y=randi([-250 250]);
        coords(:,:,end)=coords(:,:,end) + repmat([(1e-9)*pos_x (1e-9)*pos_y 0],size(coords,1),1);
        for t=time_step:time_step:ObservationTime
            m_y=randi([-max_x,max_x]);
            %m_y=randi([-max_y,max_y]);
            r=roots([1, 0, a, b - m_y^2]);
            r1 = r(imag(r)==0);
            msize = numel(r1);
            m_x=r1(randperm(msize, 1));
            if(max_z>0)
                m_z=randi([-max_z,max_z]);
            else
                m_z=0;
            end
            coords(:,:,end+1)=coords(:,:,end) + repmat([(1e-9)*m_x (1e-9)*m_y (1e-9)*m_z],size(coords,1),1);
        end
    case upper('Diffuse')
        tt=time_step:time_step:ObservationTime;
        need_check=true;
        gg=random('Normal',0,sqrt(2*DiffuseCoef*time_step),[size(coords,1),size(coords,2), length(tt)]);
        while(need_check)
            if isempty(find(coords(:,3,end)~=0))
                gg(:,3,:)=0;
            end
            pp=cumsum(gg,3);
            particles_coords_temp=repmat(coords,[1,1,length(tt)])+pp;
            xyz_xyz_span=repmat(xyz_span(:).',[size(coords,1),1,length(tt)]);
            indices_to_change=[];
            %             indices_to_change=find(abs(particles_coords_temp)>xyz_xyz_span/2);
            %             indices_to_change=find(sum(...
            %                 (abs(particles_coords_temp(:,1,:))>xyz_span(1)/2)...
            %                 | (abs(particles_coords_temp(:,2,:))>xyz_span(2)/2)...
            %                 | (abs(particles_coords_temp(:,3,:))>xyz_span(3)/2)...
            %                 ,3));
            length(indices_to_change)/length(particles_coords_temp(:))
            if isempty(indices_to_change)
                need_check=false;
            else
                %                 gg(indices_to_change,:,:) = random('Normal',0,sqrt(2*DiffuseCoef*time_step),[length(indices_to_change),size(coords,2), length(tt)]);
                gg(indices_to_change) = random('Normal',0,sqrt(2*DiffuseCoef*time_step),size(indices_to_change));
            end
        end
        coords(:,:,1+ (1:size(particles_coords_temp,3)))=particles_coords_temp;
        
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

emission=abs(round(random('Normal',photonsperstep,sqrt(photonsperstep),[NoOfEmitters NoOfTimeSteps])));
moment_bleach=NoOfTimeSteps-random('geo',p_bleach,[NoOfEmitters,1]);

% if strcmpi('Poisson_ExpBleach',Blinking.name)
%     c_bleach=repmat(exp(-(1:NoOfTimeSteps).*time_step./tau_bleach),NoOfEmitters,1);
%     c_emm=rand(size(c_bleach));
%     e_emm=c_emm>c_bleach;
%     for n=1:NoOfEmitters; 
%         a=find(e_emm(n,:),1,'first'); 
%         if isempty(a) 
%             moment_bleach(n)=NoOfEmitters+1; 
%         else
%             moment_bleach(n)=a; 
%         end
%     end
% else
    moment_bleach=NoOfTimeSteps-random('geo',p_bleach,[NoOfEmitters,1]);
% end
blinking=zeros([NoOfEmitters NoOfTimeSteps]);
moment_bleach(moment_bleach<0)=NoOfTimeSteps;
% switch upper(Blinking.name)
%     
%     case upper('Poisson')
        if flag_emission_same
            n=1;
            t_end=0;
            while t_end<NoOfTimeSteps
                T_chunk=random('Poisson',round(tau_blink/time_step));
                t_on_chunk=min(random('Poisson',T_chunk*p_on),T_chunk);
                the_on_samples=randperm(T_chunk,t_on_chunk);
                blinking(n,t_end+the_on_samples)=1;
                t_end=t_end+T_chunk;
            end
            blinking(n,moment_bleach(n):end)=0;
            blinking=repmat(blinking(n,:),NoOfEmitters,1);
            emission=repmat(emission(1,:),NoOfEmitters,1);
        else
        for n=1:NoOfEmitters
            t_end=0;
            while t_end<NoOfTimeSteps
                T_chunk=random('Poisson',round(tau_blink/time_step));
                t_on_chunk=min(random('Poisson',T_chunk*p_on),T_chunk);
                the_on_samples=randperm(T_chunk,t_on_chunk);
                blinking(n,t_end+the_on_samples)=1;
                t_end=t_end+T_chunk;
            end
            blinking(n,moment_bleach(n):end)=0;
        end
        end
%     case upper('Poisson_ExpBleach')
%         if flag_emission_same
%             n=1;
%             t_end=0;
%             while t_end<NoOfTimeSteps
%                 T_chunk=random('Poisson',round(tau_blink/time_step));
%                 t_on_chunk=min(random('Poisson',T_chunk*p_on),T_chunk);
%                 the_on_samples=randperm(T_chunk,t_on_chunk);
%                 blinking(n,t_end+the_on_samples)=1;
%                 t_end=t_end+T_chunk;
%             end
%             blinking(n,moment_bleach(n):end)=0;
%             blinking=repmat(blinking(n,:),NoOfEmitters,1);
%             emission=repmat(emission(1,:),NoOfEmitters,1);
%         else
%             for n=1:NoOfEmitters
%                 t_end=0;
%                 while t_end<NoOfTimeSteps
%                     T_chunk=random('Poisson',round(tau_blink/time_step));
%                     t_on_chunk=min(random('Poisson',T_chunk*p_on),T_chunk);
%                     the_on_samples=randperm(T_chunk,t_on_chunk);
%                     blinking(n,t_end+the_on_samples)=1;
%                     t_end=t_end+T_chunk;
%                 end
%                 blinking(n,moment_bleach(n):end)=0;
%             end
%         end
%         
%     case upper('PowerPoisson')
%         for n=1:NoOfEmitters
%             t_end=0;
%             while t_end<NoOfTimeSteps
%                 T_chunk=random('gp',round(tau_blink/time_step));
%                 t_on_chunk=min(random('Poisson',T_chunk*p_on),T_chunk);
%                 the_on_samples=randperm(T_chunk,t_on_chunk);
%                 blinking(n,t_end+the_on_samples)=1;
%                 t_end=t_end+T_chunk;
%             end
%             blinking(n,moment_bleach(n):end)=0;
%         end
% end
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