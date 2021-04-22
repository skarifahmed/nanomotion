%% Read me
% This function generates the complete 3-D PSF of the system, it takes into account the
% aberration and other factors, all. It uses dyadic Green's function and
% property of incoherence. It can also consider off-focus ccd as well as sample. 
% Other simplified versions have also been coded - see functions
% generate_PSF_<xxxx>

% Intellectual property
% This code is created by Krishna Agarwal (BioSym, SMART, MIT) on 26th Sept
% 2014. The copyright of the code is held by Krishna Agarwal and SMART. The
% code can be used freely for non-commercial research purposes. Kindly cite
% the paper below if you publish the results generated using this code and
% acknowledge the copyright holders.
%!!!!!!!!!!

% Inputs:
% NA: numerical aperture of the system
% n_obj_n_ccd: refractive indices of the focal and ccd region. For oil immersion, first value is the refractive index of oil. The second value, n_ccd is most often 1
% f_obj_f_ccd: focal lengths of the objective and ccd in microns. Generally, f_ccd/f_obj=magnification
% lam_um: wavelength of observation in micrometer
% vec_r_foc: N_foc x 3 matrix: [r theta phi] cooridnates of the sample
%           points, r is in micrometers
% vec_r_ccd: N_ccd x 3 matrix: [r theta phi] cooridnates of the ccd's pixel
%           points, r is in micrometers

% Outputs:
% I_map: intensity PSF from focus [0 0 0] to the input ccd points
% G: dyadic PSF from focus [0 0 0] to the input ccd points
% I_out: mapping from the input sample points to the input CCD points


function [I_map,G,I_out]=generate_PSF_SOFI (NA,n_obj_n_ccd,f_obj_f_ccd,lam_um,vec_r_foc,vec_r_ccd,ORDER)
%constant and non-integrands
lam=lam_um*1e-6;% changing wavelength into meters
mu0=4*pi*1e-6; % H/m permeability of free space
c=3e8;% speed of light : unit meter/s
n_obj=n_obj_n_ccd(1);
n_ccd=n_obj_n_ccd(2);
f_obj=f_obj_f_ccd(1);
f_ccd=f_obj_f_ccd(2);
theta_max_obj=asin(NA/n_obj);
vec_r_foc(:,1)=vec_r_foc(:,1)*1e-6;
vec_r_ccd(:,1)=vec_r_ccd(:,1)*1e-6;

k_ccd=2*pi*n_ccd/lam;
k_obj=2*pi*n_obj/lam;
beta=(-1i*k_ccd/(8*pi))*(f_obj/f_ccd)*sqrt(n_obj/n_ccd)*exp(i*(k_obj*f_obj + k_ccd*f_ccd));
beta2=i*2*pi*(c/lam)*mu0;

kx_ccd=k_ccd*(vec_r_ccd(:,1).*cos(vec_r_ccd(:,2)).*cos(vec_r_ccd(:,3)));
ky_ccd=k_ccd*(vec_r_ccd(:,1).*cos(vec_r_ccd(:,2)).*sin(vec_r_ccd(:,3)));
kz_ccd=k_ccd*(vec_r_ccd(:,1).*sin(vec_r_ccd(:,2)));
kx_foc=k_obj*(vec_r_foc(:,1).*cos(vec_r_foc(:,2)).*cos(vec_r_foc(:,3)));
ky_foc=k_obj*(vec_r_foc(:,1).*cos(vec_r_foc(:,2)).*sin(vec_r_foc(:,3)));
kz_foc=k_obj*(vec_r_foc(:,1).*sin(vec_r_foc(:,2)));
%integration variable and related quantities
th_obj=linspace(0,theta_max_obj,200);
th_ccd=asin((f_obj/f_ccd)*sin(th_obj));

%creating 3D matrices for different quantities_needed
[X_ccd,X_foc,Th_obj]=meshgrid(kx_ccd,kx_foc,th_obj);
[Y_ccd,Y_foc,Th_ccd]=meshgrid(ky_ccd,ky_foc,th_ccd);
[Z_ccd,Z_foc,~]=meshgrid(kz_ccd,kz_foc,th_obj);

X=-1*(X_ccd.*sin(Th_ccd) + X_foc.*sin(Th_obj)); clear X_ccd X_obj;
Y=-1*(Y_ccd.*sin(Th_ccd) + Y_foc.*sin(Th_obj)); clear Y_ccd Y_obj;
Z=Z_ccd.*cos(Th_ccd)-Z_foc.*cos(Th_obj); clear Z_ccd Z_obj;

Rho=sqrt(X.^2+Y.^2);
Psi=atan2(Y,X);
clear X Y;

A=sqrt(cos(Th_obj)./cos(Th_ccd)).*sin(Th_obj).*exp(i*Z).*th_obj(2); clear Z;

Pi_0_cc=1+cos(Th_ccd).*cos(Th_obj);
I_0_cc=sum(A.*Pi_0_cc.*besselj(0,Rho),3); clear Pi_0_cc;
Pi_0_ss=sin(Th_ccd).*sin(Th_obj);
I_0_ss=sum(A.*Pi_0_ss.*besselj(0,Rho),3); clear Pi_0_ss;
Pi_1_sc=cos(Th_ccd).*sin(Th_obj).*cos(Psi);
I_1_sc=sum(A.*Pi_1_sc.*besselj(1,Rho),3); clear Pi_1_sc;
Pi_1_ss=cos(Th_ccd).*sin(Th_obj).*sin(Psi);
I_1_ss=sum(A.*Pi_1_ss.*besselj(1,Rho),3); clear Pi_1_ss;
Pi_1_cc=sin(Th_ccd).*cos(Th_obj).*cos(Psi);
I_1_cc=sum(A.*Pi_1_cc.*besselj(1,Rho),3); clear Pi_1_cc;
Pi_1_cs=sin(Th_ccd).*cos(Th_obj).*sin(Psi);
I_1_cs=sum(A.*Pi_1_cs.*besselj(1,Rho),3); clear Pi_1_cs;
Pi_2_cc=(1-cos(Th_ccd).*cos(Th_obj)).*cos(2*Psi);
I_2_cc=sum(A.*Pi_2_cc.*besselj(2,Rho),3); clear Pi_2_cc;
Pi_2_cs=(1-cos(Th_ccd).*cos(Th_obj)).*sin(2*Psi);
I_2_cs=sum(A.*Pi_2_cs.*besselj(2,Rho),3); clear Pi_2_cs;

Gxx=beta*beta2*(I_0_cc + I_2_cc);
Gyy=beta*beta2*(I_0_cc - I_2_cc);
Gzz=beta*beta2*I_0_ss;
clear I_0_cc I_2_cc I_0_ss;
Gxy=beta*beta2*I_2_cs;
Gyx=beta*beta2*I_2_cs;
clear I_2_cs;
Gxz=beta*beta2*(-2)*i*I_1_sc;
clear I_1_sc;
Gyz=beta*beta2*(-2)*i*I_1_ss;
clear I_1_ss;
Gzx=beta*beta2*(2)*i*I_1_cc;
clear I_1_cc;
Gzy=beta*beta2*(2)*i*I_1_cs;
clear I_1_cs;

% Intensity PSF for incoherent source case:
% !!!!! check carefully
Ixx=(abs(Gxx).^2+abs(Gxy).^2+abs(Gxz).^2);
Iyy=(abs(Gyx).^2+abs(Gyy).^2+abs(Gyz).^2);
Izz=(abs(Gzx).^2+abs(Gzy).^2+abs(Gzz).^2);

G{1,1}=Gxx;clear Gxx;
G{1,2}=Gxy;clear Gxy;
G{1,3}=Gxz;clear Gxz;
G{2,1}=Gyx;clear Gyx;
G{2,2}=Gyy;clear Gyy;
G{2,3}=Gyz;clear Gyz;
G{3,1}=Gzx;clear Gzx;
G{3,2}=Gzy;clear Gzy;
G{3,3}=Gzz;clear Gzz;

I_out=(Ixx+Iyy+Izz);
I_map{1}=Ixx;clear Ixx;
I_map{2}=Iyy;clear Iyy;
I_map{3}=Izz;clear Izz;

if ~(size(vec_r_foc,1)==1)
% if ~(size(vec_r_foc,1)==1 && vec_r_foc(1)==0)
    vec_r_foc=[0.0 0 0];
[I_map,G,~]=generate_PSF_SOFI (NA,n_obj_n_ccd,f_obj_f_ccd,lam_um,vec_r_foc,vec_r_ccd*1e6,ORDER);
end