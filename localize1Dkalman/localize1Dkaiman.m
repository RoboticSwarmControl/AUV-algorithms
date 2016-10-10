%Demostration of localization 1D Kalman 1 sensor
% 
% Input: Change n to add or reduce the number of sensor 
% 
%

% Author: Haoran Zhao, Master, Electrical engineering
% University of Houston
% email address: zhaohaorandl@gmail.com 
% 
% December 1992; Last revision: 22-Sep-2016

%------------- BEGIN CODE --------------
%%------------ main variable -----------
n = 4; % sensor position
boat_init = 0;
boat_current = 0;
v= 3; %m/s
dt= .1; 
dist = zeros(n,1);
distm = zeros(n,100);
index=ones(n,1);
mean_est = zeros(n,1);
cov_m = zeros(n,1);
j_last=1;
row = zeros(n,1);
x_m = zeros(n,1);
x_new=zeros(n,1);
cov_new=zeros(n,1);
% time = zeros(1,13001);
diff = zeros(n,13001);

%%------------ random sensor position -------
a = 0;
b = 1000;
x_est = a+(b-a).*rand(n,1);
y = zeros(n,1);
sensor_position=zeros(n,1);
for j=1:n
    sensor_position(j)=x_est(j)+randn(1,1)*30;
end

dr = 200;
dr_error=zeros(n,1);
%%------------ plot sensor pdf -----------
figure(1)
subplot(2,1,1)
hold on
hBpts = plot (boat_init,'>','markerfacecolor','y','markersize',10);
grid on
hRline = plot ([boat_init(1,1)-dr boat_init(1,1)+dr],[0 0],'g');
hTtitle = title(['T = ',num2str(0)]); 
hSpts = plot (sensor_position(:,1),y(:,1),'m+');
[hMspdf, hMspts]= initmpdf(n,mean_est,cov_m);
[hKalpdf,hKalpts]=initkalpdf(n,x_new,cov_new);
[hEstpdf,hEstpts]=initestpdf(n,x_est,ones(n,1).*30);
hold off
axis([0 1000 -0.02 0.05]);

%%------------ figure(2) ------------------
dif = abs(x_est - sensor_position);
[row,~]=find(dif==max(dif));
x_m = ceil(dif(row,1))+1;
diff(:,1)=dif;
subplot(2,1,2)
hold on
[hDfpts,hDfline]=initdfdist(n,0,diff);
hold off
axis([0 400 0 x_m]);
title ('Distance bwt x-kalman and sensor-position')

%%------------ measure and postion estimate--------------
for t = 0:dt:400
    boat_current = boat_current + v*dt;
    bt = t - floor(t);
    if t==0
    col =1;
    else
    col = t/0.1+1;
    end
    col=round(col);
    diff(:,col)=dif(:,1);
    
    if bt < 10^-3
    for j = 1:n
    dist(j) = errorm([sensor_position(j),0],[boat_current,0]);
    dr_error(j) = dr+randn(1,1)*2;
    end
    
    for k = 1:n       
        if dist(k)<dr_error(k);
            d=boat_current - x_est(k);
            if d<=0
            x_m(k) = boat_current + dist(k);
            else
            x_m(k) = boat_current - dist(k);
            end
            if index(k)==1
            [x_new(k),cov_new(k)]=kalmanfilter(x_est(k),30,x_m(k));
            index(k)=0;
            else
            [x_new(k),cov_new(k)]=kalmanfilter(x_est(k),cov_new(k),x_m(k));    
            end
            x_est(k)=x_new(k);
            cov_m(k,1) = 10; 
            dif(k) = abs(x_est(k)-sensor_position(k));
            diff(k,col)=dif(k);
        else
        x_m(k)=0;
        cov_m(k,1)=0;  
        end
    end
    end
   [hDfpts,hDfline]=upddfdist(n,t,diff,hDfpts,hDfline);
   [hMspdf,hMspts]=updmpdf(n,x_m,cov_m,hMspdf,hMspts); 
   [hKalpdf,hKalpts]=updkalpdf(n,x_new,cov_new,hKalpdf,hKalpts);
%    [hEstpdf,hEstpts]=updestpdf(n,mean_new,cov_new,hEstpdf,hEstpts);
    set(hRline,'xdata',[boat_current(1,1)-dr boat_current(1,1)+dr],'ydata',[0 0]);
    set(hBpts,'xdata',boat_current(:,1),'ydata',y(:,1));
    set(hTtitle,'String',['T = ', num2str(t)]);
    drawnow;
%     pause(0.1);
end
diff(:,13001)
%------------- END OF CODE --------------
