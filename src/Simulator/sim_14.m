for set_no=1:length(Sample.Set)
    set_no
    [coords{set_no,1},emission{set_no,1},NoOfEmitters{set_no,1},ExceptionString{set_no,1}]...
        =generate_sample_and_dynamics(Sample.Set(set_no).Region,...
        Sample.Set(set_no).Pattern,...
        Sample.Set(set_no).Dynamics,...
        Sample.Set(set_no).Blinking,...
        Sample.Set(set_no).NoOfEmitters,...
        Sample.timestep,...
        Sample.ObservationTime,Sample.xyz_span,flag_emission_same);
        Sample.Set(set_no).NoOfEmittersActual=size(emission{set_no},1);
%         emission{set_no}=repmat(emission{set_no}(1,1:end),size(emission{set_no},1),1);
end

coords_all=cell2mat(coords);
emission_all=cell2mat(emission);
NoOfEmitters_all=sum(cell2mat(NoOfEmitters));
GroundTruth.e_emitters=reshape(emission_all,NoOfEmitters_all,[]);
GroundTruth.x_emitters=reshape(coords_all(:,1,:),NoOfEmitters_all,[]);
GroundTruth.y_emitters=reshape(coords_all(:,2,:),NoOfEmitters_all,[]);
GroundTruth.z_emitters=reshape(coords_all(:,3,:),NoOfEmitters_all,[]);
GroundTruth.t_emitters=0:Sample.timestep:Sample.ObservationTime;

disp('Generated sample geometry, dynamics, and emissions...');
save([save_path file_str_basic underscore 'Sample'],'Sample','GroundTruth','-v7.3');
emission_all=GroundTruth.e_emitters;
coords_all(:,1,:)=GroundTruth.x_emitters;
coords_all(:,2,:)=GroundTruth.y_emitters;
coords_all(:,3,:)=GroundTruth.z_emitters;
e_max=max(max(emission_all));
plot_str={[1 0 0],[0 1 0],[1 1 0],[0 0 1]};

% Flag=0
% if(Flag==1)
%  
%  
%  for t=1:size(coords_all,3)
%     %figure_handle=figure('Visible','on');
%     figure_handle=figure;
%     hold on;
%     rectangle('Position',[-1*Sample.xyz_span(1)/2,-1*Sample.xyz_span(1)/2,Sample.xyz_span(1),Sample.xyz_span(2)],'FaceColor',[0 0 0]);
%     No_of_previous=0;
% for set_no=1:length(Sample.Set)
%  %   plot(coords{set_no}(:,1,t),coords{set_no}(:,2,t),'Color',plot_str{set_no});
%     color_set=emission_all(No_of_previous+(1:Sample.Set(set_no).NoOfEmittersActual),t)*plot_str{set_no}/e_max;
%     for e=1:size(color_set,1)
%     plot(coords_all(No_of_previous+e,1,t),coords_all(No_of_previous+e,2,t),'Color',color_set(e,:),'Marker','+','Linewidth',1.5);
%     end
%     No_of_previous=No_of_previous+Sample.Set(set_no).NoOfEmittersActual;
% end
% 
% axis equal;axis([(Sample.xyz_span(1)/2)*[-1 1] (Sample.xyz_span(2)/2)*[-1 1]]);axis ij;
% xlabel('x (m)');
% ylabel('y (m)');
% title(['t = ' num2str(GroundTruth.t_emitters(t)*1e3) ' ms']);
% figure(gcf);
%set(gcf,'Visible','off');

%Arif Comment
% Sample_Frame(t)=getframe(figure_handle);
% pause(0.1);
% close;
% end
% writerObj = VideoWriter([file_str_basic underscore 'Sample' underscore 'Video' '.avi']);
% writerObj.FrameRate = 50;
% open(writerObj);
% for f_no=1:length(Sample_Frame)
%     writeVideo(writerObj,Sample_Frame(f_no));
% end
% close(writerObj);
% disp('Generated video of sample');
% 
% end
[ System.x_sensor,System.y_sensor,System.z_sensor ] = generate_sensor( Sensor );
System.time=Sensor.timestep*(1:Sensor.no_of_frames) + Sensor.acquisition_time_offset;

[xx,yy,zz]=meshgrid(System.x_sensor,System.y_sensor,System.z_sensor);
sz_image=size(xx);
[phi,th,r]=cart2sph(-xx(:),-yy(:),zz(:));clear xx yy zz;
vec_r_ccd=[r(:)*1e6,th(:),phi(:)]; clear r th phi;

for t=1:size(GroundTruth.x_emitters,2)
[phi,th,r]=cart2sph(GroundTruth.x_emitters(:,t),GroundTruth.y_emitters(:,t),GroundTruth.z_emitters(:,t));clear xx yy zz;
vec_r_foc=[r(:)*1e6,th(:),phi(:)]; clear r th phi;
[~,~,I_out]=generate_PSF_airy(Sensor,vec_r_foc,vec_r_ccd);
Intensity(:,t)=GroundTruth.e_emitters(:,t).'*I_out;
%disp(['set 1 ' , num2str(t)]);
end
Intensity=reshape(Intensity, length(System.x_sensor), length(System.y_sensor),[]);
disp('....');
System.Image=[];
for t=1:length(System.time)
    time_indices=((GroundTruth.t_emitters>(System.time(t)-Sensor.timestep)) & (GroundTruth.t_emitters<=System.time(t)));
    if ~isempty(time_indices)
    System.Image(:,:,t)=sum(Intensity(:,:,time_indices),3);
    else
        System.Image(:,:,t)=zeros(size(Intensity,1),size(Intensity,2));
    end
    %disp(['set 1 ' , num2str(t)]);
end
%%
System.Image=System.Image/max(System.Image(:));
disp('Generated sensor image...');
save([save_path file_str_basic underscore 'Image_ori'],'Sensor','System','-v7.3');
save([save_path file_str_basic underscore 'parameters'],'Params');
%}
load([save_path file_str_basic underscore 'Image_ori']);
%% Generate tiff movie
%%{
outputFileName=[save_path file_str_basic underscore 'Image' '.tif'];
for K=1:length(System.time)
    if K==1
        imwrite(uint16((2^16 -1)*System.Image(:, :, K)), outputFileName, 'WriteMode', 'overwrite', 'Compression','none');
    else
        imwrite(uint16((2^16 -1)*System.Image(:, :, K)), outputFileName, 'WriteMode', 'append', 'Compression','none');
    end
    
end
disp('Created image into tiff file');
%}
clear Sample GroundTruth System Sensor Intensity I_out;
