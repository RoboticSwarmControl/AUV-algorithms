function [mean,hSppts,hEstpts] = initestpts(s_number,distance_error,A_position)
mean=zeros(s_number,2);
hSppts = zeros(s_number,1);
hEstpts = zeros(s_number,1);
for i = 1:s_number
    xy=pts_in(distance_error(i,:),A_position);
    xy = empty_in(xy,distance_error(i,:),A_position);
    [row,~]=size(xy);
    mean(i,:) =roundn( sum(xy)/row,0);
%     cov_mst(:,:,w) = vpa1(cov(xy));
    hSppts(i) = plot(xy(:,1),xy(:,2),'r.');
    hEstpts(i) = plot(mean(:,1),mean(:,2),'k+'); %estimate position
end