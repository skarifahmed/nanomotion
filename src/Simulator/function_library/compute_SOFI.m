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
% Image=Image/max(Image(:));
% Image=Image*4/7;
Normalize=1;N_time;
%% Initializing Corr and intermediate variables
tau=0;-1:1;

%%
if strcmpi(SOFI.type,'correlation')
    switch SOFI.order
        case 2
            for t1=1:length(tau)
                Image1=circshift(Image,[0 tau(t1)]);
                Corr(:,t1)=mean(Image.*Image1,2)/Normalize;
            end
        case 3
            for t1 = 1:length(tau)
                Image1=circshift(Image,[0 tau(t1)]);
                for t2 = 1:length(tau)
                    Image2=circshift(Image,[0 tau(t2)]);
                    Corr(:,t1,t2)=mean(Image.*Image1.*Image2,2)/Normalize;
                end
            end
        case 4
            for t1 = 1:length(tau)
                Image1=circshift(Image,[0 tau(t1)]);
                for t2 = 1:length(tau)
                    Image2=circshift(Image,[0 tau(t2)]);
                    for t3 = 1:length(tau)
                        Image3=circshift(Image,[0 tau(t3)]);
                        Corr(:,t1,t2,t3)=mean(Image.*Image1.*Image2.*Image3,2)/Normalize;
                    end
                end
            end
        case 5
            for t1 = 1:length(tau)
                Image1=circshift(Image,[0 tau(t1)]);
                for t2 = 1:length(tau)
                    Image2=circshift(Image,[0 tau(t2)]);
                    for t3 = 1:length(tau)
                        Image3=circshift(Image,[0 tau(t3)]);
                        for t4 = 1:length(tau)
                            Image4=circshift(Image,[0 tau(t4)]);
                            Corr(:,t1,t2,t3,t4)=mean(Image.*Image1.*Image2.*Image3.*Image4,2)/Normalize;
                        end
                    end
                end
            end
        case 6
            for t1 = 1:length(tau)
                Image1=circshift(Image,[0 tau(t1)]);
                for t2 = 1:length(tau)
                    Image2=circshift(Image,[0 tau(t2)]);
                    for t3 = 1:length(tau)
                        Image3=circshift(Image,[0 tau(t3)]);
                        for t4 = 1:length(tau)
                            Image4=circshift(Image,[0 tau(t4)]);
                            for t5 = 1:length(tau)
                                Image5=circshift(Image,[0 tau(t5)]);
                                Corr(:,t1,t2,t3,t4,t5)=mean(Image.*Image1.*Image2.*Image3.*Image4.*Image5,2)/Normalize;
                            end
                        end
                    end
                end
            end
        otherwise
            for t1=1:length(tau)
                Image1=circshift(Image,[0 tau(t1)]);
                Corr(:,t1)=sum(Image.*Image1,2)/Normalize;
            end
    end
    
    Output=reshape(Corr,[sz(1:2) length(tau)*ones(1,SOFI.order-1)]);
    return;
end
%%
if strcmpi(SOFI.type,'Cumulant')
    if SOFI.order>=2
        for t1=1:length(tau)
            Image1=circshift(Image,[0 tau(t1)]);
            Corr{1}(:,t1)=mean(Image.*Image1,2)/Normalize;
        end
        Cum{1}=Corr{1};% C2=G2
    end
    if SOFI.order>=3
        for t1 = 1:length(tau)
            Image1=circshift(Image,[0 tau(t1)]);
            for t2 = 1:length(tau)
                Image2=circshift(Image,[0 tau(t2)]);
                Corr{2}(:,t1,t2)=mean(Image.*Image1.*Image2,2)/Normalize;
            end
        end
        Cum{2}=Corr{2};%C3=G3
    end
    
    if SOFI.order>=4
        for t1 = 1:length(tau)
            Image1=circshift(Image,[0 tau(t1)]);
            for t2 = 1:length(tau)
                Image2=circshift(Image,[0 tau(t2)]);
                for t3 = 1:length(tau)
                    Image3=circshift(Image,[0 tau(t3)]);
                    Corr{3}(:,t1,t2,t3)=mean(Image.*Image1.*Image2.*Image3,2)/Normalize;
                end
            end
        end
%         Corr_tau1=repmat(Corr{2},[1 1 length(tau) length(tau)]);
%         Corr_tau2=repmat(Corr{2},[1 length(tau) 1 length(tau)]);
%         Corr_tau3=repmat(Corr{2},[1 length(tau) length(tau) 1]);
%         Cum{3}=Corr{3} - Corr_tau1.*Corr_tau2 - Corr_tau2.*Corr_tau3 - Corr_tau3.*Corr_tau1;
        C2=repmat(Corr{1},[1 1 length(tau) length(tau)]);
        Cum{3}=Corr{3} - 3*C2.*C2;%Dertinger_2009_supp
%         Cum{3}=Corr{3} + 3*C2.*C2;%positive made
        %C4 = G4 - (3C2)*C2*G2;
        
    end
    if SOFI.order>=5
        for t1 = 1:length(tau)
            Image1=circshift(Image,[0 tau(t1)]);
            for t2 = 1:length(tau)
                Image2=circshift(Image,[0 tau(t2)]);
                for t3 = 1:length(tau)
                    Image3=circshift(Image,[0 tau(t3)]);
                    for t4 = 1:length(tau)
                        Image4=circshift(Image,[0 tau(t4)]);
                        Corr{4}(:,t1,t2,t3,t4)=mean(Image.*Image1.*Image2.*Image3.*Image4,2)/Normalize;
                    end
                end
            end
        end
        C3=repmat(Cum{2},[1 1 1 length(tau) length(tau)]);
        G2=repmat(Corr{1},[1 1 length(tau) length(tau) length(tau)]);
        G3=repmat(Corr{2},[1 1 1 length(tau) length(tau)]);
        C2=repmat(Cum{1},[1 1 length(tau) length(tau) length(tau)]);
        Cum{4}=Corr{4} - 6*C3.*G2  - 4*C2.*G3 ;%Dertinger_2009_supp
%         Cum{4}=Corr{4} + 6*C3.*G2  + 4*C2.*G3 ;%positive
    % C5 = G5 - (4C2)C3G2 -(4C3)C2G3
    end
    if SOFI.order>=6
        for t1 = 1:length(tau)
            Image1=circshift(Image,[0 tau(t1)]);
            for t2 = 1:length(tau)
                Image2=circshift(Image,[0 tau(t2)]);
                for t3 = 1:length(tau)
                    Image3=circshift(Image,[0 tau(t3)]);
                    for t4 = 1:length(tau)
                        Image4=circshift(Image,[0 tau(t4)]);
                        for t5 = 1:length(tau)
                            Image5=circshift(Image,[0 tau(t5)]);
                            Corr{5}(:,t1,t2,t3,t4,t5)=mean(Image.*Image1.*Image2.*Image3.*Image4.*Image5,2)/Normalize;
                        end
                    end
                end
            end
        end
        C4=repmat(Cum{3},[1 1 1 1 length(tau) length(tau)]);
        G2=repmat(Corr{1},[1 1 length(tau) length(tau) length(tau) length(tau)]);
        C3=repmat(Cum{2},[1 1 1 length(tau) length(tau) length(tau)]);
        G3=repmat(Corr{2},[1 1 1 length(tau) length(tau) length(tau)]);
        C2=repmat(Cum{1},[1 1 length(tau) length(tau) length(tau) length(tau)]);
        G4=repmat(Corr{3},[1 1 1 1 length(tau) length(tau)]);
%         Cum{5}=Corr{5} - 10*C4.*G2  - 5*C3.*G3 -10*C2.*G4 ;
        Cum{5}=Corr{5} - 10*C4.*G2  - 10*C3.*G3 -5*C2.*G4 ;%Dertinger_2009_supp
%         Cum{5}=Corr{5} + 10*C4.*G2  + 10*C3.*G3 + 5*C2.*G4 ;
    end
    Output=reshape(Cum{SOFI.order-1},[sz(1:2) length(tau)*ones(1,SOFI.order-1)]);
    return;
end
%%
%{
%Simple Second order
for n_pix=1:size(Image,1)
    tt(n_pix,:)=xcorr(Image(n_pix,:),Image(n_pix,:));
end
%Way1
figure;
subplot(1,2,1);imagesc(reshape(sum(Image,2),sz(1:end-1)));axis equal;axis tight; colormap hot;
subplot(1,2,2);imagesc(reshape(tt(:,ceil(size(tt,2)/2)),sz(1:end-1)));axis equal;axis tight; colormap hot;
%Way2
figure;
subplot(1,2,1);imagesc(reshape(sum(Image,2),sz(1:end-1)));axis equal;axis tight; colormap hot;
subplot(1,2,2);imagesc(reshape(sum(Image.^2,2),sz(1:end-1)));axis equal;axis tight; colormap hot;
Output=sum(Image.^2,2);
%}
