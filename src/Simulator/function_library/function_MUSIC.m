%% Read me
% This function generates the MUSIC pseudospectrum from a given data matrix
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
% Pseudospectrum: A (1 x M_prime array). The points with asymptotically 
%   high values represent the locations with scatterers
% M : Control parameter, used to define the kernel

%%

function [Pseudospectrum,M,signal_strength,s,d_PN,d_PS] = function_MUSIC (Data_CCD_illum, Mapping_CCD_test_pt,M,option_scale_singular,option_Pseudospectrum,window_weights,SNR_cutoff)

if isempty(window_weights)
    window_weights=eye(size(Data_CCD_illum,1));
elseif min(size(window_weights))
    window_weights=diag(window_weights);
end

Mapping_CCD_test_pt=window_weights*Mapping_CCD_test_pt/max(Mapping_CCD_test_pt(:));
Data=window_weights*Data_CCD_illum;
[u,s,v]=svd(Data,'econ');
V=abs(fftshift(fft(v),1));

s=diag(s);
if isempty(option_scale_singular)
    option_scale_singular=A;
end
switch option_scale_singular
    case 'A'
        scale_val=1;
        plot_str=[];
    case 'B'
        scale_val=s(1);
        plot_str='/ \sigma_1';
    case 'C'
        scale_val=sqrt(sum(s.^2));
        plot_str='/ ||{\bf S}||';
    case 'SNR'
        scale_val=1;
        plot_str='';
    otherwise
        scale_val=value;
        plot_str=['/' num2str(value)];
end
if isempty(M)
figure;plot(log10(s/scale_val),'*');ylabel(['log_{10}(\sigma_n' plot_str ')']); xlabel('n')
M=input('How many svds to consider as comprising of the signal subspace');
end
% figure;plot(log10(s/scale_val),'*');ylabel(['log_{10}(\sigma_n' plot_str ')']); xlabel('n')
if option_scale_singular=='SNR'
    M=find(log10(s)<SNR_cutoff,1,'first');
    if isempty(M)
        M=0;
    end
end


%{
mean_V=mean(V);
std_V=abs(std(v));

plot_V=max(V,[],1)./sum(V,1);
% for vv=1:size(V,2)
%     plot_V(vv)=sum(V(V(:,1)>mean_V(vv),vv))/sum(V(:,vv));
% end
judge_V=(cumsum(plot_V))/max(sqrt(cumsum(plot_V)));
sss=diag(s);sss(sss<1e-4)=1e-4;
% figure('units','normalized','outerposition',[0 0 1 1]);subplot(1,3,1);imagesc(V);colorbar;subplot(1,3,2);plotyy([0 1:length(s)],[0 judge_V],1:length(s),log10(sss));subplot(1,3,3);plotyy([0 1:length(s)],[0 std_V],1:length(s),log10(sss));
% plot((sum(abs(fftshift(fft(v),1)))));subplot(1,3,3);plot(log10(diag(s)),'*');

% M=find(judge_V<=0.7,1,'last');if isempty(M), M=1;end
% M=find(std_V>=0.7,1,'last');if isempty(M), M=1;end
%}

if (M >= min(size(Data))) || (M<=0)
%     if strcmpi(option_scale_singular,'c')
%         aa=find(log10(s/scale_val)<-3);
%         if isempty(aa)
%             M=min(size(Data_CCD_illum))-3;
%         else
%         M=aa(1);
%         end
%     else
%     M=min(size(Data_CCD_illum))-3;
%     end
    M=min(size(Data))-3;
end


u_red=u(:,1:M);
u_orth=u(:,(M+1):end);
s_orth=s((M+1):end);
% MM=diag(1./s_orth)*u_orth.'*Mapping_CCD_test_pt;
MM=u_orth.'*Mapping_CCD_test_pt;
NN=sqrt(sum(u_orth.^2,1)).'*sqrt(sum(Mapping_CCD_test_pt.^2,1));
PP=u_red.'*Mapping_CCD_test_pt;



switch lower(option_Pseudospectrum)
    case 'nullproj'
        Pseudospectrum=sqrt(sum((abs(MM)).^2,1)); 
    case 'sigproj'
        Pseudospectrum=sqrt(sum((abs(PP)).^2,1));
        
    case 'invnullproj'
        Pseudospectrum=1./sqrt(sum((abs(MM)).^2,1)); 
    case 'invnullprojpower2'
        Pseudospectrum=(1./sqrt(sum((abs(MM)).^2,1))).^2;
    case 'ratsigprojnullproj'
        d_s_byd_N=sqrt(sum((abs(PP)).^2,1))./sqrt(sum((abs(MM)).^2,1));
        Pseudospectrum=d_s_byd_N;
    case 'ratsigprojnullprojpower2'
        d_s_byd_N=sqrt(sum((abs(PP)).^2,1))./sqrt(sum((abs(MM)).^2,1));
        Pseudospectrum=d_s_byd_N.^2;
    case 'ratsigprojnullprojpower4'
        d_s_byd_N=sqrt(sum((abs(PP)).^2,1))./sqrt(sum((abs(MM)).^2,1));
        Pseudospectrum=d_s_byd_N.^4;
    case 'acosnullproj'
        Pseudospectrum=mean(acos(MM./NN));
    case 'asinnullproj'
        Pseudospectrum=mean(asin(MM./NN));
    case 'atannullproj'
        Pseudospectrum=mean(atan(MM./NN));
%     case 'asecnullproj'
%         Pseudospectrum=mean(acos(1./(MM./NN)));
%     case 'acosecnullproj'
%         Pseudospectrum=mean(asin(1./(MM./NN)));
    case 'acotnullproj'
        Pseudospectrum=mean(atan(1./(MM./NN)));
    case 'atanhnullproj'
        Pseudospectrum=mean(atanh(MM./NN));
    case 'invacosnullproj'
        Pseudospectrum=1./mean(acos(MM./NN));
    case 'invasinnullproj'
        Pseudospectrum=1./mean(asin(MM./NN));
    case 'invatannullproj'
        Pseudospectrum=1./mean(atan(MM./NN)); 
    case 'invacotnullproj'
        Pseudospectrum=1./mean(1./(atan(MM./NN))); 
    case 'invatanhnullproj'
        Pseudospectrum=1./mean(atanh(MM./NN));
    otherwise
        if isempty(strfind(lower(option_Pseudospectrum),'ratsigprojnullprojpower'))
        Pseudospectrum=1./sqrt(sum((abs(MM)).^2,1));
        else
            power_MUSIC=str2num(option_Pseudospectrum((length('ratsigprojnullprojpower')+1):end));
            d_s_byd_N=sqrt(sum((abs(PP)).^2,1))./sqrt(sum((abs(MM)).^2,1));
            Pseudospectrum=d_s_byd_N.^power_MUSIC;
        end
            
end
d_PN=sqrt(sum((abs(MM)).^2,1));d_PS=sqrt(sum((abs(PP)).^2,1));

% Pseudospectrum=1./sqrt(sum((abs(MM)).^2,1)); % option 1
% Pseudospectrum=sqrt(sum((abs(PP)).^2,1))./sqrt(sum((abs(MM)).^2,1)); % option 2

signal_strength=[sqrt(trace(Data_CCD_illum*Data_CCD_illum')) sqrt(sum(s(:).^2)) sqrt(sum(s(1:M).^2)) s(1)];
% pause(1);close;