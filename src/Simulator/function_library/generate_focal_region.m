%% Read me
% This function generates the grid for focal region and the distribution
% that determines the probability of finding an emitter

% Intellectual property
% This code is created by Krishna Agarwal (BioSym, SMART, MIT) on 25th Sept
% 2014. The copyright of the code is held by Krishna Agarwal and SMART. The
% code can be used freely for non-commercial research purposes. Kindly cite
% the following paper if you publish the results generated using this code and
% acknowledge the copyright holders.
% !!!!!

% Inputs:
% FocalRegion: a structure containing the following elements:
%   FocalRegion.lateral_view_area : view area's lateral cross-section is a square and this number gives teh side of the square in miicrometers
%   FocalRegion.lateral_element_size : for computation purposes, the view area should also be discretized. 'lateral_element_size' indicates the size of each grid element in micrometers along x- and y- axes
%   FocalRegion.longitudinal_view_area : view area's longitudinal extent is given be this variable in miicrometers
%   FocalRegion.longitudinal_element_size : for computation purposes, the view area should also be discretized. 'longitudinal_element_size' indicates the size of each grid element along the z-axis in micrometers
%   FocalRegion.distribution_map : parameters to generate the probability distribution map - ie the probability at which emitters may be present here.
%       FocalRegion.distribution_map.name : type of distribution function
%           (options and parameters given below)
%       FocalRegion.distribution_map.param : cells - each cell containing
%           one parameter of the distribution function
%       =====================================================
%       function            Parameters
%       =====================================================
%       constant/default    []
%       gaussian            mu (3-elements vector), sig (3-elements vector)
%       sphere              r: (micrometer)
%       rbf1 [type: 1./(r.^n + a.^n)] a: (micrometer), n: order
%       rbf2 [type: 1./((r+r1).^n + (r-r2).^n + a.^n)]
%                           r1,r2,a: (micrometer), n:order
%       ellipsoid           a,b,c: micrometer
%       shell               r_out,r_in: (micrometer)
%       shell1              r_out,r_in: (micrometer), p: p1/p2, where p1 is
%                           probability of finding emitter in the inner
%                           layer and p2 is the probability of finding
%                           emitter in the outer layer
%       NPointCircle        use only for 2D, r (micrometer): radius of the circle, N :
%                           number of emitters uniformly distributed
%       NPointMCircle       use only for 2D, r (micrometer): radius of each circle, N :
%                           number of emitters uniformly distributed in
%                           each circle, R: radius of the bigger circle over
%                           which smaller scircles are distributed, M: no
%                           of circles
%       NPointLine          l (micrometer): total extent of the x-directed line, N :
%                           number of emitters uniformly distributed
%       NPointMConceptricCircle   use only for 2D. {R1,N1,R2,N2,R3,N3...}
%       NPointGrid          a (micrometer): size of square/cube, N :
%                           number of emitters uniformly distributed along
%                           each direction
%       NPointHex           use only for 2D, a (micrometer): size of one edge of hexagon, N :
%                           number of emitters uniformly distributed along
%                           each side of hexagon (N>=2)
%       NPointAll           All points as emitter locations - useful for
%                           generating PSF for MUSIC testing
%   odd_or_even: 'odd' or 'even' or [], input which specifies if odd-numbered
%       grid or even-numbered grid is generated

% Outputs:
% x_foc,y_foc,z_foc : 1D vectors containing the x,y,z-coordinates of the
%                       focal region's computation grid.
% N_x_foc,N_y_foc,N_z_foc: number of grid points along x,y,z-directions
% distribution_map: probability distribution map of the presence of
%                   emitter.

%%
function [x_foc,y_foc,z_foc,N_x_foc,N_y_foc,N_z_foc,distribution_map]=generate_focal_region(FocalRegion,odd_or_even,no_of_pixels)

N_x_foc = ceil(FocalRegion.lateral_view_area/FocalRegion.lateral_element_size);
%generates even numbered focal grid
if mod(N_x_foc,2)==1 & strcmpi(odd_or_even,'even')
    N_x_foc=N_x_foc+1;
end
%generates odd numbered pixel grid
if mod(N_x_foc,2)==0 & strcmpi(odd_or_even,'odd')
    N_x_foc=N_x_foc+1;
end
if ~isempty(no_of_pixels)
    N_x_foc=no_of_pixels;
end
    
x_foc = linspace(-0.5*FocalRegion.lateral_view_area,0.5*FocalRegion.lateral_view_area,N_x_foc+1);
x_foc = (x_foc(1:end-1) + x_foc(2:end))/2;
N_y_foc = N_x_foc;
y_foc = x_foc;
N_z_foc = ceil(FocalRegion.longitudinal_view_area/FocalRegion.longitudinal_element_size);
z_foc = linspace(-0.5*FocalRegion.longitudinal_view_area,0.5*FocalRegion.longitudinal_view_area,N_z_foc+1);
z_foc = (z_foc(1:end-1) + z_foc(2:end))/2;

distribution_map=generate_distribution_map(FocalRegion.distribution_map,x_foc,y_foc,z_foc);

end

%%
function distribution_map = generate_distribution_map(Distribution_parameters,x_foc,y_foc,z_foc)
% [xx,yy,zz]=meshgrid(round(x_foc/(x_foc(2)-x_foc(1))),round(y_foc/(y_foc(2)-y_foc(1))),round(z_foc/(z_foc(2)-z_foc(1))));
% the above tries to assign integers to grid - may not work well if the
% element sizes are different in different directions
if length(z_foc)>1
    [xx,yy,zz]=meshgrid(x_foc/(x_foc(2)-x_foc(1)),y_foc/(y_foc(2)-y_foc(1)),z_foc/(z_foc(2)-z_foc(1)));
    xx=round(xx);
    yy=round(yy);
    zz=round(zz);
else
    [xx,yy,zz]=meshgrid(x_foc/(x_foc(2)-x_foc(1)),y_foc/(y_foc(2)-y_foc(1)),z_foc);
    xx=round(xx);
    yy=round(yy);
end
distribution_map = zeros(size(xx));
switch lower(Distribution_parameters.name)
    case 'constant'
        distribution_map = ones(size(xx));
    case 'gaussian'
        if (length(z_foc)>1)
            mu=Distribution_parameters.param{1}./[(x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)),(z_foc(2)-z_foc(1))];
            sig=Distribution_parameters.param{2}./[(x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)),(z_foc(2)-z_foc(1))].^2;
        else
            mu=Distribution_parameters.param{1}./[(x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)),1];
            sig=Distribution_parameters.param{2}./[(x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)),1].^2;
        end
        %         sig=Distribution_parameters.param{2};
        distribution_map = mvnpdf([xx(:) yy(:) zz(:)],mu,diag(sig));
        distribution_map = reshape(distribution_map,size(xx));
    case 'sphere'
        r = Distribution_parameters.param{1}/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));
        distribution_map(xx.^2+yy.^2+zz.^2<r.^2)=1;
    case 'rbf1'
        a=Distribution_parameters.param{1}/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));
        n=Distribution_parameters.param{2};
        r=sqrt(xx.^2+yy.^2+zz.^2);
        distribution_map=1./(r.^n + a.^n);
    case 'rbf2'
        a=Distribution_parameters.param{1}/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));
        n=Distribution_parameters.param{2};
        r1=Distribution_parameters.param{3}/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));
        r2=Distribution_parameters.param{4}/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));
        r=sqrt(xx.^2+yy.^2+zz.^2);
        distribution_map=abs(1./((r+r1).^n + (r-r2).^n + a.^n));
    case 'ellipsoid'
        a=Distribution_parameters.param{1}/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));
        b=Distribution_parameters.param{2}/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));
        c=Distribution_parameters.param{3}/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));
        distribution_map((xx./a).^2+(yy./b).^2+(zz./c).^2<1)=1;
    case 'shell'
        r_out = max(Distribution_parameters.param{1},Distribution_parameters.param{2})/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));;
        r_in = min(Distribution_parameters.param{1},Distribution_parameters.param{2})/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));;
        distribution_map((xx.^2+yy.^2+zz.^2<r_out.^2) & (xx.^2+yy.^2+zz.^2>r_in.^2))=1;
    case 'shell1'
        r_out = max(Distribution_parameters.param{1},Distribution_parameters.param{2})/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));;
        r_in = min(Distribution_parameters.param{1},Distribution_parameters.param{2})/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));;
        distribution_map((xx.^2+yy.^2+zz.^2<r_out.^2) & (xx.^2+yy.^2+zz.^2>r_in.^2))=1;
        distribution_map((xx.^2+yy.^2+zz.^2<min(r_out,r_in).^2) )= Distribution_parameters.param{3} ;
    case 'npointcircle'
        N=Distribution_parameters.param{2};
        th=(2*pi/N)*(1:N);
        r=Distribution_parameters.param{1}/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));
        x=round(r*cos(th));
        y=round(r*sin(th));
        xx=round(xx);
        yy=round(yy);
        for ii=1:N
            dist=sqrt((xx-x(ii)).^2 + (yy-y(ii)).^2);
            [~,ind]=min(dist(:));
            distribution_map(ind)=1;
            %             distribution_map(xx==x(ii) & yy==y(ii))=1;
        end
    case 'npointmcircle'
        N=Distribution_parameters.param{2};
        M=Distribution_parameters.param{4};
        r=Distribution_parameters.param{1}/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));
        R=Distribution_parameters.param{3}/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));
        th=(2*pi/N)*(1:N);
        Th=(2*pi/M)*(1:M);
        x0=R*cos(Th);
        y0=R*sin(Th);
        
        xx=round(xx);
        yy=round(yy);
        for mm=1:M
            x=round(r*cos(th)-x0(mm));
            y=round(r*sin(th)-y0(mm));
            for ii=1:N
                dist=sqrt((xx-x(ii)).^2 + (yy-y(ii)).^2);
                [~,ind]=min(dist(:));
                distribution_map(ind)=1;
                %             distribution_map(xx==x(ii) & yy==y(ii))=1;
            end
        end
        
    case 'npointline'
        N=Distribution_parameters.param{2};
        l=Distribution_parameters.param{1}/(x_foc(2)-x_foc(1));
        x=round(linspace(-l/2,l/2,N));
        for ii=1:N
            dist=sqrt((xx-x(ii)).^2 + yy.^2);
            [~,ind]=min(dist(:));
            distribution_map(ind)=1;
            %             distribution_map(xx==x(ii) & yy==min(abs(yy(:))))=1;
        end
        
    case 'npointliney'
        N=Distribution_parameters.param{2};
        l=Distribution_parameters.param{1}/(y_foc(2)-y_foc(1));
        y=round(linspace(-l/2,l/2,N));
        for ii=1:N
            dist=sqrt(xx.^2 + (yy-y(ii)).^2);
            [~,ind]=min(dist(:));
            distribution_map(ind)=1;
            %             distribution_map(xx==x(ii) & yy==min(abs(yy(:))))=1;
        end
        
    case 'npointmliney'
        M=Distribution_parameters.param{end};
        sep_lines=linspace(-Distribution_parameters.param{end-1}/2,Distribution_parameters.param{end-1}/2,M)/(y_foc(2)-y_foc(1));
        
        DD=reshape(Distribution_parameters.param,2,[]);
        
        for jj=1:M
            N=DD{2,jj};
            l=DD{1,jj}/(y_foc(2)-y_foc(1));
            y=round(linspace(-l/2,l/2,N));
        for ii=1:N
            dist=sqrt((xx-sep_lines(jj)).^2 + (yy-y(ii)).^2);
            [~,ind]=min(dist(:));
            distribution_map(ind)=1;
            %             distribution_map(xx==x(ii) & yy==min(abs(yy(:))))=1;
        end
        end
        
    case 'npointgrid'
        N=Distribution_parameters.param{2};
        l=Distribution_parameters.param{1};
        x=round(linspace(-l/2,l/2,N)/(x_foc(2)-x_foc(1)));
        y=round(linspace(-l/2,l/2,N)/(y_foc(2)-y_foc(1)));
        for ii=1:N
            for jj=1:N
                dist=sqrt((xx-x(ii)).^2 + (yy-y(jj)).^2);
                [~,ind]=min(dist(:));
                distribution_map(ind)=1;
                %             distribution_map(xx==x(ii) & yy==y(jj))=1;
            end
        end
    case 'npointmconcentriccircle'
        for circle_no=1:(length(Distribution_parameters.param)/2);
            r=Distribution_parameters.param{(2*circle_no)-1}/max((x_foc(2)-x_foc(1)),(y_foc(2)-y_foc(1)));
            N=Distribution_parameters.param{(2*circle_no)};
            th=(2*pi/N)*(1:N);
            x=round(r*cos(th));
            y=round(r*sin(th));
            
            xx=round(xx);
            yy=round(yy);
            for ii=1:N
                dist=sqrt((xx-x(ii)).^2 + (yy-y(ii)).^2);
                [~,ind]=min(dist(:));
                distribution_map(ind)=1;
            end
        end
    case 'npointhex'
        distribution_map = ones(size(xx));
        % will develop later
    case 'npointall'
        distribution_map = ones(size(xx));
        % will develop later
    otherwise
        distribution_map = ones(size(xx));
end
distribution_map=distribution_map/max(distribution_map(:));
distribution_map=distribution_map/sum(distribution_map(:));
end