function distances=cursor_distances(cc,dimensions)

for ii=1:length(cc)
    coords(ii,:)=cc(ii).Position;
end

for dim=1:dimensions
    [x1,x2]=meshgrid(coords(:,dim),coords(:,dim))
    dist_2(:,:,dim)=(x1-x2).^2;
end
distances=sqrt(sum(dist_2,3));