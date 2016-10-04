%motion planning to aid localize underwater sensor in minmum time and energy cost
%
%
%Underwter communication have to be considered bunch of error during
%process. Here we only consider the motion planning of the boat, so we
%don't precisely compute the error in the process. Thus, assume error
%concludes random error and system error, which can be consider as Gaussian
%distribution N(1, 0.033) for anchor sensor 
% 
%The algorithm is based on range based schemes localization
%
%
% Outline:   
%    random generate n fixed sensor's position
%    approximately compute initial position of each sensor by anchor sensor
%    use minimal spanning tree algorithm (MST) 
%    use delaunay triangulation algorithm (DT) 
%    use voronoi disgrams algorithm (VD)
%    add adjusted kalman filter to each algorithm 
%    compare the time cost vs error ratio
%    compare the path length vs error ratio
%    compare the energy cost of boat and acoustic communication vs error ratio
%    compare the final area of covariance in same time, path longth, energy cost of each algorithm  
%
% result: 
%    
%
%
%
% 
% Subfunctions: dist_error.m dist.m G_error.m initplotellipse.m
% initPlotmvn.m overlapping.m pts_in.m updplotellipse.m updplotmvn.m vpa1.m
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
s_number = 1 ; % number of sensor
n=s_number;
v = 1; % speed of boat
boat_init = [20,30];
dt=1; %sampling rate
boat_last_mst = boat_init(:,:);
boat_new_mst = boat_init(:,:);
boat_current_mst = boat_init(:,:);
distance_mst = zeros(s_number,1);
A_position = [0 0;0 100; 100 100; 100 0]; %archor beacon for initial locolization
distance_error=zeros(s_number,4);
r = 400; % sensor detect range
Ad =1; % anchor odometry error sigma
Sd=0.5; %sensor odometry error sigma
% r= r+randn*5;
index=1;
recal = 1;
bsdist=zeros(s_number,1);
bsdist_error=zeros(s_number,1);
%-------------- Random generate sensor position (Actual position) -------------
% a = 0;
% b = 100;
% c = 0;
% d = 100;
% x = a+(b-a).*rand(n,1);
% y = c+(d-c).*rand(n,1);
% sensor_position = [x,y];
 sensor_position = [63,76]; %27,44;

%-------------- Init estimate position of each sensor ---------------
for q = 1:s_number
    distance_error(q,:)=roundn(A_error(sensor_position(q,:),A_position,Ad,Sd),-1);
end
%%---------- Plot background in figure(1) --------------
figure(1)
xlabel('x'), ylabel('y');
axis([0 100 0 100]);
hTtitle1 = title(['T = ',num2str(0)]); 
grid on 
hold on
[mean_mst,hSppts,hEstpts] = initestpts(s_number,distance_error,A_position); %%Equal posibility sample point
hSspts = plot(sensor_position(:,1),sensor_position(:,2),'g+');
hBtpts = plot(boat_current_mst(1,1), boat_current_mst(1,2),'>','markerfacecolor','y','markersize',5); %boat current position
hTrojac = plot([boat_last_mst(1,1),boat_current_mst(1,1)],[boat_last_mst(1,2), boat_current_mst(1,2)],'g'); % boat trjactory
hold off

%%---------- build world ---------------------------
xrange = 0:100; %initialize world (allspace is 0)
yrange = 0:100;
[X,Y]=ndgrid(xrange,yrange);
range = [X(:),Y(:)];
belief= zeros(101,101);
for i = 1:s_number
    xy{i}=pts_in(distance_error(i,:),A_position);
    xy{i} = empty_in(xy{i},distance_error(i,:),A_position);
    particlemap=xy{i};
    %initialize probabilities
    [numParticle,~]=size(xy{i});
    probs{i} = ones(numParticle,1)/numParticle; % uniformly distributed
    for i1=1:numParticle
    belief(particlemap(i1,1),particlemap(i1,2))=probs{i}(i1);
    end
end


%%---------- plot background in figure(2) -------------
figure(2)
hTtitle2 = title(['T = ',num2str(0)]);
redrawWorld(belief);

%%---------- main process -----------------
for t= 1:dt:120
    %--------------- Recursive localiation based on MST -----------------
     if recal ==1
         for y=1:s_number
            distance_mst(y,1)=dist(boat_last_mst(1,:),mean_mst(y,:));
         end
         distance_mst1=distance_mst;
         distance_mst1=sort(distance_mst1);
         recal =0;
         status =1;
     end
   %%-------------- compute nearest sensor -------------- 
    if status ==1
        [row1,~]=find(distance_mst==distance_mst1(index,1));
        boat_new_mst = mean_mst(row1,:);
%%------------- compute boat bearing -----------------        
        theta = atan2(mean_mst(1,2)-boat_last_mst(1,2),mean_mst(1,1)-boat_last_mst(1,1));
        v_mst = [v*cos(theta),v*sin(theta)];
    end
    
%%------------- compute boat current position --------
boat_current_mst = [boat_current_mst(1,1)+v_mst(1,1)*dt,boat_current_mst(1,2)+v_mst(1,2)*dt];
%%----- compute distance bwt boat and sensor position ----------
delta = dist(boat_current_mst,boat_new_mst);
if delta < v*dt
    v_new = delta/dt;
    v_new_mst = [v_new*cos(theta),v_new*sin(theta)];
    boat_current_mst = boat_current_mst+v_new_mst.*dt;
    boat_last_mst = boat_new_mst;
end
% if delta < v*dt
%     v_new = delta/dt;
%     v_new_mst = [v_new*cos(theta),v_new*sin(theta)];
%     boat_current_mst = boat_current_mst+v_new_mst.*dt;
%     boat_last_mst = boat_new_mst;
%      %%----------- update drawing figure(1) --------------
%     set(hTrojac, 'xdata',[boat_current_mst(1,1), boat_last_mst(1,1)], 'ydata',[boat_current_mst(1,2), boat_last_mst(1,2)]);
%     set(hBtpts, 'xdata', boat_current_mst(1,1), 'ydata', boat_current_mst(1,2));
%     set(hTtitle1,'String',['T = ', num2str(t)]);
% %     set(hTtitle2,'String',['T = ', num2str(t)]);
%     drawnow;
% else
    



for j= 1:s_number
bsdist(j) = dist(sensor_position(j,:),boat_current_mst);
bsdist_error(j) = bsdist+randn(1,1)*Sd;
if bsdist < 50
    particlemap=xy{j};
    %initialize probabilities
    [numParticle,~]=size(particlemap);
    pbdist=zeros(numParticle,1);
    % compute dist from boat to particle
    for i1=1:numParticle
    pbdist(i1) = dist(boat_current_mst,xy{i}(i1,:));
    end
    
    % update probility of each particle
    for i2=1:numParticle
        probs{j}(i2)=updtprob(pbdist(i2),bsdist_error(j),Sd);    
    end
    
    normalizer = 1/sum(sum(probs{j}));
    probs{j}=probs{j}.*normalizer;
    
    for i3=1:numParticle
    belief(particlemap(i3,1),particlemap(i3,2))=probs{j}(i3);
    end
    
    
end
end


%%----------- update drawing figure(1) -------------
for i4=1:s_number
    m=probs{i4};
    n=xy{i4};
    fidx=m>10^-3;
    idx=find(fidx);
    [ridx,~]=size(idx);
    xy_new=zeros(ridx,2);
    for i5=1:ridx
        xy_new(i5,:)=n(idx(i5),:);
    end
    set(hSppts(i4),'xdata',xy_new(:,1),'ydata',xy_new(:,2));
end
    set(hTrojac, 'xdata',[boat_current_mst(1,1), boat_last_mst(1,1)], 'ydata',[boat_current_mst(1,2), boat_last_mst(1,2)]);
    set(hBtpts, 'xdata',boat_current_mst(1,1),'ydata',boat_current_mst(1,2));
    set(hTtitle1,'String',['T = ', num2str(t)]);
%     set(hTtitle2,'String',['T = ', num2str(t)]);
    drawnow;

%%--------- update figure(2)-------------
    figure(2)
    redrawWorld(belief);
    pause(.3);
end






%------------- END OF CODE --------------