%% Read me
% This function generates the svd from a given data matrix
% Data_CCD_illum and a given mapping matrix Mapping_CCD_test_pt

% Intellectual property
% This code is created by Krishna Agarwal (BioSym, SMART, MIT) on 21st Nov
% 2014. The copyright of the code is held by Krishna Agarwal and SMART. The
% code can be used freely for non-commercial research purposes. Kindly cite
% the paper below if you publish the results generated using this code and
% acknowledge the copyright holders.

%% References of relevance:
% Agarwal, Krishna, and Xudong Chen. "Applicability of MUSIC-type imaging 
% in two-dimensional electromagnetic inverse problems." Antennas and 
% Propagation, IEEE Transactions on 56.10 (2008): 3217-3223.

% Chen, Xudong, and Krishna Agarwal. "MUSIC algorithm for two-dimensional 
% inverse problems with special characteristics of cylinders." Antennas and
% Propagation, IEEE Transactions on 56.6 (2008): 1808-1812.

% Agarwal, Krishna, et al. "Practical applications of multiple signal 
% classification." International Journal of RF and Microwave Computer?Aided 
% Engineering 22.3 (2012): 359-369.

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
% s : singular values

%%

function [s] = function_SVD (Data_CCD_illum,window_weights)

if isempty(window_weights)
    window_weights=eye(size(Data_CCD_illum,1));
elseif min(size(window_weights))
    window_weights=diag(window_weights);
end

Data=window_weights*Data_CCD_illum;
s=svd(Data,'econ');
