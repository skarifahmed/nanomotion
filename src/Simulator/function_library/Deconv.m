clc;clear all; close all;
D_ph=[50 25 10];
% D_ph=[50 10];
PSF_zeros=[0.19 0.16 0.15 0.15]
PSF_FWHMs=[0.1549 0.1179 0.1077 0.1067]
lam=1.340;%wavelngth in micron

% TL_TR_BL_x=[
%     295 590 390 485
%     645 940 740 830
%     300 590 390 475
%     ];
% TL_TR_BL_y=[
%     400 495 475 275
%     395 500 475 275
%     745 845 815 620
%     ];

% TL_TR_BL_x=[
%     290 590 390 480
%     640 940 740 825
%     295 590 390 475
%     ];
% TL_TR_BL_y=[
%     400 495 475 275
%     395 495 475 275
%     745 845 815 620
%     ];

%240 nm pitch, 50 nm pitch not mentioned
TL_TR_BL_x=[
    290 590 390 480
    635 940 735 825
    295 590 390 475
    ];
TL_TR_BL_y=[
    400 495 475 275
    395 495 475 275
    740 845 810 620
    ];

%240 nm pitch, 50 nm pitch mentioned
TL_TR_BL_x=[
    281 590 390 480
    627 940 735 825
    280 590 390 475
    ];
TL_TR_BL_y=[
    379 495 475 275
    381 495 475 275
    728 845 810 620
    ];


%%% 400 nm ptch
% TL_TR_BL_x=[
%     212 272
%     789 852
%     216 269
%     ];
% TL_TR_BL_y=[
%     109 398
%     109 400
%     681 975
%     ];


base_points=[
50 50
950 50
50 950
];

indices_to_slice=[500:900];[400:450 850:950];

view_area=82;%microns
zoom=4;
no_of_pixels = 1024;


pixel_size=view_area./(zoom*no_of_pixels);
PSF_zero = PSF_zeros*lam;
N=ceil(PSF_zero*2./pixel_size);
PSF_FWHM = PSF_FWHMs*lam;
S=(PSF_FWHM./(2*pixel_size));

for ii=1:length(D_ph)
    PSF = fspecial('gaussian',N(ii),S(ii));PSF=PSF/max(PSF(:));
%     PSF = ones(N(ii));PSF=PSF/max(PSF(:));
    figure;
%     subplot(1,3,1);imagesc(pixel_size*linspace(-1,1,N(ii)),pixel_size*linspace(-1,1,N(ii)),PSF);caxis([0 max(PSF(:))]);axis square;
    im=(imread(['Im' num2str(D_ph(ii)) '.tif']));
    [J P] = deconvblind(im,PSF);
    J=double(J);
    J=J/max(J(:));
    im=double(im);
    im=im/max(im(:));
    subplot(2,2,1);imshow(im);subplot(2,2,2);imshow(J);title(['Pinhole Diameter ' num2str(D_ph(ii))]);
    tform=cp2tform([TL_TR_BL_x(:,ii) TL_TR_BL_y(:,ii)], (base_points), 'affine');
    im_f=imtransform(im,tform,'Xdata',[1 1000],'Ydata',[1 1000]);
    J_f=imtransform(J,tform,'Xdata',[1 1000],'Ydata',[1 1000]);
    J_f=J_f/max(J_f(:));
    im_f=im_f/max(im_f(:));
    subplot(2,2,3);imshow(im_f);subplot(2,2,4);imshow(J_f);title(['Pinhole Diameter ' num2str(D_ph(ii))]);
    Set1{ii}=(J_f(indices_to_slice,1:500)).';
    Set2{ii}=J_f(1:500,indices_to_slice);
end
%%
figure;hold on;
str={'r','g','b','c'}
for ii=1:length(D_ph)
    plot_data=1*mean(Set1{ii},2);
    plot_data=(plot_data-min(plot_data));
    plot_data=plot_data/max(plot_data);
    plot(plot_data,str{ii});
end
figure;hold on;
for ii=1:length(D_ph)
    plot_data=1*mean(Set2{ii},2);
    plot_data=(plot_data-min(plot_data));
    plot_data=plot_data/max(plot_data);
    plot(plot_data,str{ii});
end