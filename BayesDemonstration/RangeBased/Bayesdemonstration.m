%function BayesDemonstration()
%BayesDemonstration is used to demonstrate Bayes' rule for an underwater
%localization sestup
%
%The algorithm uses range-based localization
%
%
% Subfunctions: updatprob, drawworld, redrawworld, drawcontour
%
% Author: Haoran Zhao, Master, Electrical Engineering
% University of Houston
% email address: zhaohaorandl@gmail.com
% Last revision: 4-Oct-2016

%------------- BEGIN CODE --------------
clear
clc
digits(4);
format compact
%%------------ main variables -----------
sensor_position = [0,0];
boat_position=zeros(1,2);
xrange=-50:50; %define an area
yrange=-50:50;
[X,Y]=ndgrid(xrange,yrange);
xy=[X(:),Y(:)];
r=numel(X);
sd=1;   %standard deviation of range measurements
md=1;
step = 5;
range = 30;
dist_cell=zeros(r,1);
probs = zeros(r,1);
belief_init = ones(101,101)./numel(ones(101,101)); % normalize initial probability
belief=belief_init;
pts_in = cell(1,step);
dn_belief = cell(1,step);
b_position=zeros(101,101);
b_var=NaN(1,step);
%%----------- main code ---------------


for i=1:step
    %%% random generate boat position
    theta = rand(1,1)*2*pi;
    sbdist = (range+md)* rand(1,1);
    boat_position = [sbdist*cos(theta),sbdist*sin(theta)];
    dist_error = sbdist+randn(1,1);
    % Update probability distribution
    for i1=1:r
        dist_cell(i1,1)=dist(xy(i1,:),boat_position); % measure distance bwt boat and sensor
        probs(i1)=updtprob(dist_cell(i1),dist_error,sd);
    end
    % Map donut sample point probability to area
    for i2=1:r
        dn_belief{i}(xy(i2,2)+51,xy(i2,1)+51)=probs(i2);
    end
    % Normalize donut probability
    normalizerb=1/sum(dn_belief{i}(:));
    Pdonut=dn_belief{i}.*normalizerb; % probability of each time estimated donut
    b_position(round(boat_position(1,2))+51,round(boat_position(1,1))+51)=max(Pdonut(:));
    % Apply the numerator of Bayes rule
    %p(a_t|b_t)=n*p(b_t|a_t-1)*p(a_t-1)
    for i4=1:i
        belief=belief.*dn_belief{i};
    end
    % Apply the denominator of Bayes rule
    %%%  m=m.*belief_init.^(i-1);
    normalizer = 1/sum(belief(:));
    belief=belief.*normalizer;
    
    %%% find possible postion in measured distance
    fidx1=probs>10^-5;
    idx1=find(fidx1);
    [r1,~]=size(idx1);
    for i3=1:r1
        pts_in{i}(i3,:) = xy(idx1(i3,1),:) ;
    end
    
    %%% compute variance of estimate position
    fidx2=belief>10^-5;
    [idr,idc]=find(fidx2==1);
    b_var(1,i)=var(idr,idc);
    %%% -----------Plotting code------------
    
    %     figure(1)
    subplot(2,2,1)
    hPpts = plot(pts_in{i}(:,1),pts_in{i}(:,2),'r.');
    hold on
    axis equal
    hSpts = plot(sensor_position(1,1),sensor_position(1,2),'b+');
    hBpts = plot(boat_position(1,1),boat_position(1,2),'^','markerfacecolor','y','Markersize',10);
    hold off
    axis([-50 50 -50 50]);
    %     legend('estimate position','sensor position','boat position');
    title({['Histogram after ',num2str(i),'  distance measurements'],['Histogram variance = ', num2str(b_var(1,i))]})
    
    %     figure(2)
    subplot(2,2,2)
    drawcontour(belief);
    title({['Histogram after ',num2str(i),'  distance measurements'],['Histogram variance = ', num2str(b_var(1,i))]})
    
    % %     figure(3)
    %     subplot(2,2,3)
    %     redrawWorld(belief);
    %     hold on
    %     redrawWorld(Pdonut);
    %     redrawWorld(b_position);
    %     title({['Histogram after ',num2str(i),'  distance measurements'],['Histogram variance = ', num2str(b_var)]})
    %
    %     figure(4)
    subplot(2,2,3)
    redrawWorlds(belief);
    hold on
    redrawWorlds(Pdonut);
    redrawWorlds(b_position);
    title({['Histogram after ',num2str(i),'  distance measurements'],['Histogram variance = ', num2str(b_var(1,i))]})
    
    subplot(2,2,4)
    C=hsv(1);
    plot(i,b_var(1,i),'.','color',C(1,:));
    hold on
    plot(1:i,b_var(1,1:i),'color',C(1,:));
    axis([1 step 0 300])
    title(['Variance after ',num2str(i),'th  range detection'])
    
    pause(1);
end
%%% display final estimate result
[r_est,c_est]=find(belief==max(belief(:)));
r_est=r_est-51;
c_est=c_est-51;
display(['Estimate of position: (x,y)= [',num2str(r_est),',',num2str(c_est),']',', step =',num2str(i)])

%------------- END OF CODE --------------