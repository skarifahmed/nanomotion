%Inputs:
% SOFI: structure containing two elements
%   SOFI.type: We can choose from the following options
%       1. Cumulant
%       2. Correlation
%       3. Balanced
%   SOFI.order: scalar, order of cumulant or correlation
%Image : the image stack captured by the imaging system
%       4-D, N_x_ccd x N_y_ccd x N_z_ccd x N_t_ccd

function Output = compute_SOFI(SOFI,Image)

% correlation functions have to be computed irrespective of the type of SOFI
% !!!! Decision: what tau to use?

%% Preprocessing
sz=size(Image);
N_time=sz(end);
Image_0=reshape(Image,[],N_time);
Image=Image_0 - repmat(mean(Image_0,2),1,N_time);
Image=Image/max(Image(:));
Normalize=N_time;
%% Initializing Corr and intermediate variables
tau=-1:1;

%%
% Let us do it for tau=0 first;
for n = 1 : SOFI.order
    interim_str='Image';
    if n>1
    for k = 1:n
        interim_str=[interim_str '.*Image'];
    end
    end
    string_eval=['mean(' interim_str ',2);'];
    Corr{n}=eval(string_eval);
    
    Cum{n}=Corr{n};
    if n>1
    for k=1:n-1
        disp([num2str(min(Cum{n})) '.......' num2str(max(Cum{n})) '........' num2str(factorial(n-1)/(factorial(k)*factorial(n-1-k))) '.......' num2str(min(Cum{n-k}.*Corr{k}))  '.......' num2str(max(Cum{n-k}.*Corr{k}))])
        Cum{n}=Cum{n}-(factorial(n-1)/(factorial(k)*factorial(n-1-k)))*Cum{n-k}.*Corr{k};
        
    end
    end
    disp(['..................................................................................'])
    disp([num2str(n) '.......' num2str(min(Corr{n})) '.......' num2str(min(Cum{n})) '.......' num2str(max(Corr{n})) '.......' num2str(max(Cum{n})) ])
    disp(['..................................................................................'])
end

if strcmpi(SOFI.type,'correlation')
Output=reshape(Corr{n},sz(1:2));
elseif strcmpi(SOFI.type,'Cumulant')
    if SOFI.order==5 || SOFI.order==6
        Output=-1*reshape(Cum{SOFI.order},sz(1:2));
    else
        Output=reshape(Cum{SOFI.order},sz(1:2));
    end
return;
end

