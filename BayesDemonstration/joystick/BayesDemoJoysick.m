%function BayesRangefreeJoystick()
%this function is used to demonstrate Bayes's rule for range free
%underwater localization with joystick input in 2D-space. In order to select boat
%position to reduce cost and time
%
% Requires a joystick
% Author: Haoran Zhao, Master, Electrical Engineering
% Advisor & Supervisor : Dr. Becker, Assosiate Professor, Electrical Engineering
% University of Houston
% email address: zhaohaorandl@gmail.com
% Last revision: 13-Oct-2016


% NOTE: To stop the program, hit ctrl+c.
% Otherwise, it will run indefinitely.
%------------- BEGIN CODE --------------
clear
close all
format compact

%------------- Main Variables -----------
% Define joystick ID
ID = 1;
% Create joystick variable
joy=vrjoystick(ID);
s_number = 5;
b_number = 1;
%generate sensor position
a=0;
b=200;
c=b/10;
d=9*b/10;
% x1=c+(d-c)*rand(s_number,1);
% y1=c+(d-c)*rand(s_number,1);
% sensor_position=[x1(:),y1(:)];
sensor_position=[37.2874 132.1360;93.5801 159.5577;92.1412 28.3507;108.181825 55.1490;148.8647 93.5427];

%generate map area
xrange = 0:b;
yrange = 0:b;
[X,Y]=ndgrid(xrange,yrange);
xy=[X(:),Y(:)];

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
v=5; %every inout from joystick move 1, 1m/s
C = hsv(s_number);
% weight = cell(b_number,s_number);

% loop until ctrl+c
while 1
    X=axis(joy, 1);     % X-axis is joystick axis 1
    Y=axis(joy, 2);     % Y-axis is joystick axis 2
    if abs(X)>0.3
        x=boat_last(1,1)+round(X)*v;% Correlate axis to angular orientation of joystick
    else
        x=boat_last(1,1);
    end
    if abs(Y)>0.3
        y=boat_last(1,2)-round(Y)*v;      % See above.
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
                pts_in{i,i1}=overlapping(xy,boat_current(i,:),range+3*sd);
                boat_belief{i,i1}=zeros(b+1,b+1);
                if abs(X)>0.3 || abs(Y)>0.3
                    for i5 = 1: numel(pts_in{i,i1})/2
                        boat_belief{i,i1}(pts_in{i,i1}(i5,2)+1,pts_in{i,i1}(i5,1)+1)=PDF(pts_in{i,i1}(i5,:),boat_current(i,:),range,sd); %prob_in{i3,i4}
                    end
                end
                belief_last{1,i1}=belief{1,i1};
                belief{1,i1}=belief{1,i1}.*boat_belief{i,i1};
                normalizer=1/sum(belief{1,i1}(:));
                belief{1,i1}=belief{1,i1}.*normalizer;
            else
                idx_in(i,i1)=0;
            end
        end
    end
    
    [r1,c1] = find(idx_in ==1);
    
    for i=1:s_number
        [r2,c2]=find(belief{1,i}>0);
        variance(1,i)=var(r2(:),c2(:));
    end
    
    [~,c3]=find(variance==0);
    [~,c4]=find(isnan(variance));
    if isempty(c3)==0
        belief{1,c3}=belief_last{1,c3};
    elseif isempty(c4)==0
        for i=1:numel(c4)
            belief{1,c4(i)}=belief_last{1,c4(i)};
        end
    end
    
    for  ii=1:s_number
        sumr=0;
        sumc=0;
        [r4,c4]=find(belief{1,ii}>belief_init);
        for i=1:numel(r4)
            sumr=sumr+belief{1,ii}(r4(i),c4(i))*(r4(i)-1);
            sumc=sumc+belief{1,ii}(r4(i),c4(i))*(c4(i)-1);
        end
        est_position(ii,2)=sumr;
        est_position(ii,1)=sumc;
        %                 est_position(ii,2)=mean(r4)-1;
        %                 est_position(ii,1)=mean(c4)-1;
        %         est_position(ii,:)=sum(weight{1,ii});
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
    %     hFig = figure(1);
    %     set(gcf,'PaperPositionMode', 'auto')
    %     set(hFig,'Position',[100 200 1000 700])
    subplot(2,2,1)
    for i=1:b_number
        for ii=1:s_number
            if isempty(find(pts_in{i,ii})>0)==0
                plot(pts_in{i,ii}(:,1),pts_in{i,ii}(:,2),'r.');
                hold on
            end
        end
    end
    plot(boat_current(1,1),boat_current(1,2),'>','MarkerSize',10,'MarkerFaceColor','y');
    hold on
    if X_b == 1
        plot(sensor_position(:,1),sensor_position(:,2),'b+');
    end
    if B_b ==1
        plot(trajactory(:,1),trajactory(:,2),'b-');
    end
    hold off
    grid on
    %     axis  equal
    axis([0 b 0 b])
    title(['Run distance = ', num2str(distance)])
    
    %     figure(2)
    subplot(2,2,2)
    plot(boat_current(:,1),boat_current(:,2),'>','markerfacecolor','y','markersize',10);
    hold on
    for i = 1: s_number
        plot(est_position(i,1),est_position(i,2),'+','color',C(i,:));
    end
    if X_b == 1
        plot(sensor_position(:,1),sensor_position(:,2),'b+');
    end
    for ii=1:s_number
        if belief{1,ii}(1,1)~=belief_init(1,1)
            contour(belief{1,ii},5);
        end
    end
    hold off
    %     axis equal
    axis([0 b 0 b])
    title(['Coutour plot after',num2str(index-1),'th detection']);
    
    %     figure(3)
    subplot(2,2,3)
    for i7 = 1:s_number
        if isempty(find(belief{1,i7})>0)==0
            redrawWorlds(belief{1,i7});
        end
        hold on
    end
    hold off
    zlim([0 1])
    axis tight
    title(['Probability Histogram after ',num2str(index-1),'th detection']);
    
    %     figure(4)
    subplot(2,2,4)
    for i2=1:s_number
        etapts = plot(index,eta(i2,index),'.','color',C(i2,:));
        hold on
        etaline = plot (1:index,eta(i2,1:index),'color',C(i2,:));
    end
    axis([0 400 0 1])
    title(['Confidence Level after ',num2str(index-1),'th detection'])
    % delay between plots
    %     pause(.05)
    drawnow;
    
end