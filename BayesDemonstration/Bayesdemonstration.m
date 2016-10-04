%This code is use for demonstrating Bayes rule can implement to underwater
%localization
%
%
%
% 
%The algorithm is based on range based schemes localization
%
%
% Outline:   
%    set a fix sensor position
%    set five different position of boat 40 away from sensor 
%    each position measure distance bwt sensor
%    implement bayes rule
%    draw 2D image
%    draw histogram
%
%
% result: 
%    
%
%
%
% 
% Subfunctions: 
% 
%

% Author: Haoran Zhao, Master, Electrical Engineering
% University of Houston
% email address: zhaohaorandl@gmail.com
% Last revision: 26-Sep-2016

%------------- BEGIN CODE --------------
clear
clc
digits(4);
%%------------ main variable -----------
sensor_position = [0,0];
boat_position=zeros(1,2);
xrange=-50:50;
yrange=-50:50;
[X,Y]=ndgrid(xrange,yrange);
xy=[X(:),Y(:)];
[r,~]=size(xy);
r=round(r);
sd=1;
md=1;
step = 5;
range = 30;
dist_cell=zeros(r,1);
probs = zeros(r,1);
belief_init = ones(101,101)./numel(ones(101,101));
m=belief_init;
%%----------- main code ---------------
for i=1:step
    theta = rand(1,1)*2*pi;
    sbdist = (range+md)* rand(1,1);
    boat_position = [sbdist*cos(theta),sbdist*sin(theta)];
    dist_error = sbdist+randn(1,1);
    for i1=1:r
         dist_cell(i1,1)=dist(xy(i1,:),boat_position);
         probs(i1)=updtprob(dist_cell(i1),dist_error,sd);
    end
    
    for i2=1:r
        belief{i}(xy(i2,1)+51,xy(i2,2)+51)=probs(i2);
    end
    normalizerb=1/sum(belief{i}(:));
    n=belief{i}.*normalizerb;
    for i4=1:i
    m=m.*belief{i};
    end
%     m=m.*belief_init.^(i-1);
    normalizer = 1/sum(m(:));
    m=m.*normalizer;
    b_var=var(m(:));
    fidx=probs>10^-5;
    idx=find(fidx);
    [r1,~]=size(idx);
    for i3=1:r1
       pts_in{i}(i3,:) =xy(idx(i3,1),:) ;
    end
    %%----------- figure(1) ---------------
figure(1)
hPpts = plot(pts_in{i}(:,1),pts_in{i}(:,2),'r.');
hold on
axis equal
grid on
hSpts = plot(sensor_position(1,1),sensor_position(1,2),'b+');
hBpts = plot(boat_position(1,1),boat_position(1,2),'k+');
axis([-50 50 -50 50]);
legend('estimate position','sensor position','boat position');
title(['Histogram after ',num2str(i),'  distance measurements'])

figure(2)
redrawWorld(m);
hold on
redrawWorld(n);
title({['Histogram after ',num2str(i),'  distance measurements'],['Histogram variance = ', num2str(b_var)]})

figure(3)
redrawWorlds(m);
hold on
redrawWorlds(n);

title({['Histogram after ',num2str(i),'  distance measurements'],['Histogram variance = ', num2str(b_var)]})


end

[r_est,c_est]=find(m==max(m(:)));
r_est=r_est-51
c_est=c_est-51





%------------- END OF CODE --------------