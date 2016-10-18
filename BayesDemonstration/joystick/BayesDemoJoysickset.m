%function BayesRangefreeJoystick()
%this function is used to demonstrate Bayes's rule for range free
%underwater localization with joystick input in 2D-space. In order to select boat
%position to reduce cost and time
%
% Requires a joystick
% Author: Haoran Zhao, Master, Electrical Engineering
% University of Houston
% email address: zhaohaorandl@gmail.com
% Last revision: 13-Oct-2016


% NOTE: To stop the program, hit ctrl+c.
% Otherwise, it will run indefinitely.
%------------- BEGIN CODE --------------
clear
clc
format compact

%------------- Main Variable -----------
% Define joystick ID
ID = 1;
% Create joystick variable
joy=vrjoystick(ID);
s_number = 5;
b_number = 1;
%generate sensor position
a=0;
b=200;
C=b/10;
d=9*b/10;
x1=C+(d-C)*rand(s_number,1);
y1=C+(d-C)*rand(s_number,1);
sensor_position=[x1(:),y1(:)];

%generate map area
xrange = 0:b;
yrange = 0:b;
[X,Y]=ndgrid(xrange,yrange);
xy=[X(:),Y(:)];
pts_in = cell(b_number,s_number); % find the possible position in the range

range = 30; %boat sensor detect range
sd = 1; %sensor error
ttgrid = 3405;
belief = cell(1,s_number); % Bayes' rule belief
belief_init = ones(b+1,b+1)./numel(ones(b+1,b+1));

for i=1:s_number
    belief{1,i}=belief_init;
end

trajactory = NaN(10000,2);
trajactory(1,:)=[0,0];
index =1 ;
boat_belief = cell(b_number,s_number); % Boat knows if there is sensor in the range
b_position = zeros(b_number,s_number); % Boat history position
prob_in = cell (b_number,s_number);
idx_in = zeros(b_number,s_number); % record if sensor in range with 0 or 1

distance_bws=zeros(b_number,s_number);
variance = zeros(1,s_number);
est_position = zeros(s_number,2);
eta = NaN(s_number,1000); % confidence level
eta_b=NaN(b_number,s_number); %confidence level each boat and each sensor

boat_last = [0,0];
boat_current = [0,0];
distance = 0;
v=1; %every inout from joystick move 1, 1m/s
X_b=0;
B_b=0;
%------------ plot world ---------
figure(1)
subplot(2,2,1)
hPtsin = initptsin(s_number,b_number,pts_in);
hBpts = plot(boat_current(1,1),boat_current(1,2),'>','MarkerSize',10,'MarkerFaceColor','y');
hold on
grid on
if X_b == 1
    hSpts = plot(sensor_position(:,1),sensor_position(:,2),'b+');
end
if B_b ==1
    hTline = plot(trajactory(:,1),trajactory(:,2),'b-');
end
hold off
axis equal
axis([0 b 0 b])
hTtitle1 = title(['Run distance = ', num2str(distance)]);

subplot(2,2,2)
axis equal
axis([0 b 0 b])
title('contour plot');
hBpts;
hold on
hEpts = plot(est_position(:,2),est_position(:,1),'g+');
if X_b ==1
    hSpts;
end
hCline = initcontour(s_number,belief);
hold off

subplot(2,2,3)
hHistogram = inithist(s_number,belief);
hold off
zlim([0 1])
axis tight
title('Probability Histogram');

subplot(2,2,4)
[confPts, confLine] = initconf(s_number,index,eta);



%------------ Main Code ----------
% loop until ctrl+c
while 1
    X=axis(joy, 1);     % X-axis is joystick axis 1
    Y=axis(joy, 2);     % Y-axis is joystick axis 2
    if abs(X)>0.5
        x=boat_last(1,1)+round(X);% Correlate axis to angular orientation of joystick
    else
        x=boat_last(1,1);
    end
    if abs(Y)>0.5
        y=boat_last(1,2)-round(Y);      % See above.
    else
        y=boat_last(1,2);
    end
    
    if x<0
        x=0;
    elseif x>b
        x=b;
    end
    
    if y<0
        y=0;
    elseif y>b
        y=b;
    end
    
    X_b = button(joy,1); % joystick X button to plot sensor position
    A_b = button(joy,2); % joystick A button to change speed
    B_b = button(joy,3); % joystick B button to plot trajactory
    Y_b = button(joy,4); % joystick Y button to change detect range and detect error
    
    if A_b ==1
        v = input('Input speed: ');
    end
    
    if Y_b ==1
        range = input('Input detect range: ');
        sd = input('Input detect error: ');
        %-----calculate ttgrid --------
        tt=overlapping(xy,[b/2,b/2],range+3*sd);
        ttgrid = numel(tt)/2;
    end
    
    boat_current = [x,y];
    distance = distance + dist(boat_last,boat_current);
    trajactory(index+1,:)=boat_current;
    if abs(x-boat_last(1,1))>0 || abs(y-boat_last(1,2))>0
        index = index+1;
    end
    boat_last = boat_current;
    %----------- sensor detect ------------
    pts_in = cell(b_number,s_number); % find the possible position in the range
    belief_last=cell(1,s_number);
    for i=1:b_number
        for i1=1:s_number
            % calculate distace bwt boat and each sensor
            distance_bws(i,i1)=dist(sensor_position(i1,:),boat_current(i,:));
            detectrange = range+randn(1,1)*sd;
            if distance_bws(i,i1) < detectrange || distance_bws(i,i1) == detectrange
                idx_in(i,i1)=1;
            else
                idx_in(i,i1)=0;
            end
        end
    end
    
    [r1,c1] = find(idx_in ==1);
    if isempty(r1)==0
        for i2=1:numel(r1)
            %find possible position in the range
            pts_in{r1(i2),c1(i2)}=overlapping(xy,boat_current(r1(i2),:),range+3*sd);
        end
        for i3 = 1:b_number
            for i4=1:s_number
                if isempty(pts_in{i3,i4})==0
                    prob_in{i3,i4}=2/numel(pts_in{i3,i4});
                    boat_belief{i3,i4}=zeros(b+1,b+1);
                    for i5 = 1: numel(pts_in{i3,i4})/2
                        boat_belief{i3,i4}(pts_in{i3,i4}(i5,2)+1,pts_in{i3,i4}(i5,1)+1)=prob_in{i3,i4};
                    end
                    belief_last{1,i4}=belief{1,i4};
                    belief{1,i4}=belief{1,i4}.*boat_belief{i3,i4};
                    normalizer=1/sum(belief{1,i4}(:));
                    belief{1,i4}=belief{1,i4}.*normalizer;
                end
            end
        end
    end
    
    for i=1:s_number
        [r2,c2]=find(belief{1,i}>0);
        variance(1,i)=var(r2(:),c2(:));
    end
    
    [~,c3]=find(variance==0);
    [~,c4]=find(isnan(variance));
    if isempty(c3)==0
        belief{1,c3}=belief_last{1,c3};
    elseif isempty(c4)==0
        belief{1,c4}=belief_last{1,c4};
    end
    
    for  ii=1:s_number
        [r4,c4]=find(belief{1,ii}>0);
        est_position(ii,1)=mean(r4)-1;
        est_position(ii,2)=mean(c4)-1;
    end
    
    %--------------- Confidence level ---------
    if isempty(r1)==0
        for ii1=1:s_number
            if numel(find(belief{1,ii1}>belief_init))==0
                eta(ii1,index)=0;
            else
                eta(ii1,index)=1-numel(find(belief{1,ii1}>belief_init))/ttgrid;
            end
        end
    else
        if index==1
            eta(:,index)=0;
        else
            eta(:,index)=eta(:,index-1);
        end
    end
    
    %----------- Plot boat curretn position ------------------
    set(hTtitle1,'String',['Run distance = ', num2str(distance)]);
    hPtsin = updtptsin(s_number,b_number,pts_in,hPtsin);
    set(hBpts,'xdata',boat_current(1,1),'ydata',boat_current(1,2));
    if X_b ==1
        set(hSpts,'xdata',sensor_position(:,1),'ydata',sensor_position(:,2));
    end
    if B_b ==1
        set(hTline,'xdata',trajactory(:,1),'ydata',trajactory(:,2));
    end
    set(hEpts,'xdata',est_position(:,2),'ydata',est_position(:,1));
    hCline = updtcontour(s_number,belief,hCline);
    hHistogram = updthist(s_number,belief,hHistogram);
    [confPts, confLine] = updtconf(s_number,index,eta,confPts,confLine);
    drawnow;
    
    % delay between plots
    drawnow
end