%%

% Rows: CCD pixels
% columns: appear to be X blocks - X=length(x_foc)
% Each block has Y columns - Y=length(y_foc)

% Assmues the folowing form of PSF coordinates:
% [
% [(r1_ccd,x1y1),(r1_ccd,x1y2),...,(r1_ccd,x1yY)],
% [(r1_ccd,x2y1),(r1_ccd,x2y2),...,(r1_ccd,x2yY)],
% ...
% [(r1_ccd,xXy1),(r1_ccd,xXy2),...,(r1_ccd,xXyY)];
% [(r2_ccd,x1y1),(r2_ccd,x1y2),...,(r2_ccd,x1yY)],
% [(r2_ccd,x2y1),(r2_ccd,x2y2),...,(r2_ccd,x2yY)],
% ...
% [(r2_ccd,xXy1),(r2_ccd,xXy2),...,(r2_ccd,xXyY)];
% and so on


%%
function [G_PSF,x_foc1,y_foc1]=interp_G_PSF(PSF_file_name,scale_interp)



load(PSF_file_name)
step_size=1/scale_interp;
X=length(x_foc);
Y=length(y_foc);
x_foc1=interp1(1:length(x_foc),x_foc,1:step_size:X);
y_foc1=interp1(1:length(y_foc),y_foc,1:step_size:Y);

X_new=length(x_foc1);
Y_new=length(y_foc1);

% expnasion in each block
for xx=1:X
    for CCD=1:size(G_PSF,1)
        G_PSF1{CCD,xx}=interp1(1:Y,G_PSF(CCD,(xx-1)*X + (1:Y)),1:step_size:Y);
    end
end
G_PSF2=cell2mat(G_PSF1);
G_PSF=G_PSF2;

clear G_PSF1 G_PSF2

% expansion between blocks
for xx=1:Y_new
    for CCD=1:size(G_PSF,1)
        G_PSF1(CCD,xx+(0:Y_new:(Y_new*(X_new-1))))=interp1(1:X,G_PSF(CCD,xx+(0:Y_new:(Y_new*(X-1)))),1:step_size:X);
end
end

G_PSF=G_PSF1;
end
