%% Read me
% This function generates the LSM pseudospectrum from a given data matrix
% Data_CCD_illum and a given mapping matrix Mapping_CCD_test_pt

% Intellectual property
% This code is created by Krishna Agarwal (BioSym, SMART, MIT) on 28th May 
% 2014. The copyright of the code is held by Krishna Agarwal and SMART. The
% code can be used freely for non-commercial research purposes. Kindly cite
% the paper below if you publish the results generated using this code and
% acknowledge the copyright holders.

%% References of relevance:


%% Inputs 
% Data_CCD_illum: Data matrix of size (N x T)
%   N: number of CCD elements 
%   T: no of independent illuminations
%   elements in Data_CCD_illum may be complex for a phase measurement
%       system OR they may be real numbers for an intensity measurement system
%   one column of Data_CCD_illum contains one CCD image (laid in raster
%       fashion)

% Mapping_CCD_test_pt: Mapping matrix (N x M_prime)
%   N: number of CCD elements
%   M_prime: number of test points
%   'm'th column in Mapping_CCD_test_pt contains the mapping from 'm'th 
%       test points to all the points on the CCD (laid in raster fashion) 
%   IF Data_CCD_illum is complex-valued, Mapping_CCD_test_pt has to be 
%       complex-valued and if Data_CCD_illum is real-valued,
%       Mapping_CCD_test_pt has to be real valued

% M: number of singular values after which to consider the noise subspace
%   if M=[], the user is asked for the value;
%   else input M is used

% option_scale_singular values: 
%     'A': no scaling
%     'B': scaling by the highest value
%     'C': scaling by net amplitude , i.e., sqrt(sum(s.^2))
%     value: scaling by value
    
    

%% Output
% Pseudospectrum: A (1 x M_prime array). The points with asymptotically 
%   high values represent the locations with scatterers
% M : Control parameter, used to define the kernel

%%

function [Pseudospectrum,U] = function_MUSIC (Data_CCD_illum, Mapping_CCD_test_pt,window_weights)
if isempty(window_weights)
    window_weights=eye(size(Data_CCD_illum,1));
elseif min(size(window_weights))
    window_weights=diag(window_weights);
end

% window_weights=window_weights.*window_weights;

Mapping_CCD_test_pt=window_weights*Mapping_CCD_test_pt/max(Mapping_CCD_test_pt(:));
Data=window_weights*Data_CCD_illum;
U=pinv(Data)*Mapping_CCD_test_pt;% minimum norm
% U=Data\Mapping_CCD_test_pt;% min non-zero elements
% sum((Mapping_CCD_test_pt-Data*U).^2,2);
%{

[u,s,v]=svd(Data); S=diag(s);
alpha=0.02;0.001;0.01;0.1;0.1;0.01;0.005;
S_inv=S./(S.^2+alpha^2);
U=v(:,1:length(S_inv))*diag(S_inv)*u'*Mapping_CCD_test_pt;
%}

%{
% TWIST
addpath('./external_function_library/TwIST_v2')
% [U,x_debias,objective,times,debias_start,mses,max_svd] = ...
%          TwIST(Mapping_CCD_test_pt,Data,0.1,'verbose',0);
for ii=1:size(Mapping_CCD_test_pt,2)
     [U(:,ii),x_debias,objective,times,debias_start,mses,max_svd] = ...
         TwIST(Mapping_CCD_test_pt(:,ii),Data,1,'verbose',0,'PHI',@(x) norm(x,1));
end
%}

%{
for ii=1:size(Mapping_CCD_test_pt,2)
    U(:,ii)=Data\Mapping_CCD_test_pt(:,ii);
end
%}

%% If Gen. Tick regularization is used
% for test_points=1:size(Mapping_CCD_test_pt,2)
%     U=pinv(Data)*
% end
%%
power=4;
Pseudospectrum_U=sum(abs(U).^power,1);sqrt(sum(U.^2,1));

Pseudospectrum_R=(sum((Mapping_CCD_test_pt-Data*U).^2,1)./(sum((Mapping_CCD_test_pt).^2,1)));
R_min=0.5;R_max=1.0;
RP_min=min(Pseudospectrum_R(:));
RP_max=max(Pseudospectrum_R(:));
Pseudospectrum_R = R_min+(R_max-R_min)*(Pseudospectrum_R-RP_min)/(RP_max-RP_min);
Pseudospectrum_R = Pseudospectrum_R.^power;

Pseudospectrum=Pseudospectrum_U.*Pseudospectrum_R;
% Pseudospectrum=Pseudospectrum_U;
% Pseudospectrum=Pseudospectrum_R;
% figure(gcf);subplot(1,2,1);imagesc(U);colorbar;subplot(1,2,2);imagesc(abs(U));colorbar;
