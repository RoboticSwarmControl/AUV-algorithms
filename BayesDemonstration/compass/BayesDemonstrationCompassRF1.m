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
boat_pass = ones(b+1,b+1); %.*100./numel(ones(b+1,b+1))
for i=1:s_number
    belief{1,i}=belief_init;
end

trajactory = NaN(10000,2);
trajactory(1,:)=[0,0];
index =1 ;
t=1;
boat_belief = cell(b_number,s_number); % Boat knows if there is sensor in the range
b_position = zeros(b_number,s_number); % Boat history position
prob_in = cell (b_number,s_number);
idx_in = zeros(b_number,s_number); % record if sensor in range with 0 or 1

distance_bws=zeros(b_number,s_number);
distance_bwest=zeros(b_number,s_number);

variance = zeros(1,s_number);
variance_last = zeros(1,s_number);
est_position = zeros(s_number,2);
eta = NaN(s_number,1000); % confidence level
eta_b=NaN(b_number,s_number); %confidence level each boat and each sensor

boat_last = [0,0];
boat_current = [0,0];
distance = 0;
v=5; %every inout from joystick move 1, 1m/s
C = hsv(s_number);
% weight = cell(b_number,s_number);
auto = 0;
idx_deltal = [];
dpmt = 0.3;
refine_m=0;
enter_refine_mode=0;
recal =1;
idx_in_last=zeros(b_number,s_number);
rsearch_mode=0;
position_r=NaN(1,2);
cl_level=0.95;
now_min=[];
now_min_old=0;
% loop until ctrl+c
while auto==0
    X=axis(joy, 1);     % X-axis is joystick axis 1
    Y=axis(joy, 2);     % Y-axis is joystick axis 2
    X_b = button(joy,1); % joystick X button to plot sensor position
    A_b = button(joy,2); % joystick A button to change speed
    B_b = button(joy,3); % joystick B button to plot trajactory
    Y_b = button(joy,4); % joystick Y button to change detect range and detect error
    if Y_b ==1
        if auto == 0
            auto =1;
        elseif auto ==1
            auto =0;
        end
        %         range = input('Input detect range: ');
        %         sd = input('Input detect error: ');
        %         %-----calculate ttgrid --------
        %         tt=overlapping(xy,[b/2,b/2],range+3*sd);
        %         ttgrid = numel(tt)/2;
    end
    
    if A_b ==1
        v = input('Input speed: ');
    end
    
    if auto ==0
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
        
        boat_current = [x,y];
        
        distance = distance + dist(boat_last,boat_current);
        trajactory(t+1,:)=boat_current;
        if abs(x-boat_last(1,1))>0 || abs(y-boat_last(1,2))>0
            t=t+1;
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
                            boat_pass(pts_in{i,i1}(i5,2)+1,pts_in{i,i1}(i5,1)+1)=0;
                        end
                    end
                    belief_last{1,i1}=belief{1,i1};
                    belief{1,i1}=belief{1,i1}.*boat_belief{i,i1};
                    normalizer=1/sum(belief{1,i1}(:));
                    belief{1,i1}=belief{1,i1}.*normalizer;
                else
                    idx_in(i,i1)=0;
                    pts_in{i,i1}=overlapping(xy,boat_current(i,:),range+3*sd);
                    if abs(X)>0.3 || abs(Y)>0.3
                        for i5 = 1: numel(pts_in{i,i1})/2
                            boat_pass(pts_in{i,i1}(i5,2)+1,pts_in{i,i1}(i5,1)+1)=0;
                        end
                    end
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
            index = index+1;
            for ii1=1:s_number
                if numel(find(belief{1,ii1}>belief_init))==0
                    eta(ii1,t)=0;
                else
                    eta(ii1,t)=1-numel(find(belief{1,ii1}>belief_init))/ttgrid;
                end
            end
        else
            if t==1
                eta(:,t)=0;
            else
                eta(:,t)=eta(:,t-1);
            end
        end
        
        %----------- Plot boat curretn position ------------------
        %     hFig = figure(1);
        %     set(gcf,'PaperPositionMode', 'auto')
        %     set(hFig,'Position',[100 200 1000 700])
        figure(1)
        subplot(2,2,1)
        colormap(jet)
        for i=1:b_number
            for ii=1:s_number
                if idx_in(i,ii)==1
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
        [anglemap,compass] = drawCompass(belief,boat_pass);
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
        drawPass(boat_pass);
        hold off
        zlim([0 1])
        axis tight
        title(['Probability Histogram after ',num2str(index-1),'th detection']);
        
        %     figure(4)
        subplot(2,2,4)
        for i2=1:s_number
            etapts = plot(t,eta(i2,t),'.','color',C(i2,:));
            hold on
            etaline = plot (1:t,eta(i2,1:t),'color',C(i2,:));
        end
        axis([0 2000/v 0 1])
        title(['Confidence Level after ',num2str(index-1),'th detection'])
        %     delay between plots
        
        pause(.001)
    end
end
%%------------------ auto trajactory planning -------------------
%%---------------------- press button Y -------------------------
while auto ==1
    if recal ==1
        [rr,cc]=find(compass==max(compass(:)));
        hp=compass(1,cc(1));
    end
    
    for t=1:1001
        X=axis(joy, 1);     % X-axis is joystick axis 1
        Y=axis(joy, 2);     % Y-axis is joystick axis 2
        X_b = button(joy,1); % joystick X button to plot sensor position
        A_b = button(joy,2); % joystick A button to change speed
        B_b = button(joy,3); % joystick B button to plot trajactory
        Y_b = button(joy,4); % joystick Y button to change detect range and detect error
        if Y_b ==1
            if auto == 0
                auto =1;
            elseif auto ==1
                auto =0;
            end
            %         range = input('Input detect range: ');
            %         sd = input('Input detect error: ');
            %         %-----calculate ttgrid --------
            %         tt=overlapping(xy,[b/2,b/2],range+3*sd);
            %         ttgrid = numel(tt)/2;
        end
        
        if A_b ==1
            v = input('Input speed: ');
        end
        
        if compass(1,cc(1))<dpmt*hp && enter_refine_mode==0
            recal = 1;
        end
        
        
        
        if enter_refine_mode ==0 && recal ==1
            [rr,cc]=find(compass==max(compass(:)));
            hp=compass(1,cc(1));
            recal = 0;
            [rr1,cc1]=find(anglemap==cc(1));
            goal = [mean(rr1)-1,mean(cc1)-1];
            lamda = atan2(goal(1,2)-boat_last(1,2),goal(1,1)-boat_last(1,1)); % heading
        end
        
        if enter_refine_mode ==1 && first ==1
            if isempty(now_min)==1
                distance_bwest=zeros(1,s_number);
                for i = 1: numel(c_r)
                    distance_bwest(1,c_r(i))=dist(boat_last,est_position(c_r(i),:));
                end
                [~,c_m]=find(distance_bwest>0);
                [~,c_m]=find(distance_bwest==min(distance_bwest(1,c_m(:))));
                gama = atan2(boat_last(1,2)-est_position(c_m,2),boat_last(1,1)-est_position(c_m,1));
                position_r(1,1)=est_position(c_m,1)+cos(gama)*range;
                position_r(1,2)=est_position(c_m,2)+sin(gama)*range;
                lamda = atan2(position_r(1,2)-boat_last(1,2),position_r(1,1)-boat_last(1,1));
                first=0;
            end
        end
        
        if enter_refine_mode == 1 && first ==1
            if isempty(now_min)==0
                gama = atan2(boat_last(1,2)-est_position(now_min,2),boat_last(1,1)-est_position(now_min,1));
                position_r(1,1)=est_position(now_min,1)+cos(gama)*range;
                position_r(1,2)=est_position(now_min,2)+sin(gama)*range;
                lamda = atan2(position_r(1,2)-boat_last(1,2),position_r(1,1)-boat_last(1,1));
                first=0;
                first_enter=1;
                c_m=now_min;
            end
        end
        
        if rsearch_mode==1 && enter_refine_mode==1
            if first_enter==1
                lamda1=atan2(boat_last(1,2)-est_position(c_m,2),boat_last(1,1)-est_position(c_m,1));
                first_enter=0;
            end
            
            x=est_position(c_m,1)+cos(lamda1-start*v/range)*range;  % boat x coordinate.
            
            y=est_position(c_m,2)+sin(lamda1-start*v/range)*range;  % boat y coordinate.
            start=start+1;
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
            
            if dist([x,y],boat_last)>v
                gama = atan2(boat_last(1,2)-est_position(c_m,2),boat_last(1,1)-est_position(c_m,1));
                position_r(1,1)=est_position(c_m,1)+cos(gama)*range;
                position_r(1,2)=est_position(c_m,2)+sin(gama)*range;
                lamda = atan2(position_r(1,2)-boat_last(1,2),position_r(1,1)-boat_last(1,1));
                rsearch_mode=0;
            end
            now_est=find(eta(:,t-1)>0);
            if eta(c_m,t-1)>cl_level
                if isempty(now_est)==0
                    now_min=find(eta(:,t-1)==min(eta(now_est(:),t-1)));
                    if now_min ~= now_min_old
                        first=1;
                    end
                    now_min_old=now_min;
                end
            end    
                if isempty(now_est)==1 || numel(now_est)==numel(find(eta(:,t-1)>cl_level))
                    enter_refine_mode=0;
                    rsearch_mode=0;
                    eta(:,t)=eta(:,t-1);
                    trajactory(t+1,:)=trajactory(t,:);
                    recal=1;
                    first=1;
                    now_min=[];
                    position_r=NaN(1,2);
                    continue
                end
                if numel(find(eta(:,t)>cl_level))==s_number  %%mean(eta(:,now))>0.95
                    break
                end
            
        end
        
        if rsearch_mode==0
            dist_bpr=dist(boat_last,position_r); %distance between boat and position_r
            x=boat_last(1,1)+cos(lamda)*v;  % boat x coordinate.
            y=boat_last(1,2)+sin(lamda)*v;  % boat y coordinate.
            if dist_bpr<v && isnan(dist_bpr)==0
                x=boat_last(1,1)+cos(lamda)*dist_bpr;
                y=boat_last(1,2)+sin(lamda)*dist_bpr;
                rsearch_mode=1;
                first_enter=1;
            end
            if x<0
                x=0;
                first=1;
            elseif x>b
                x=b;
                first=1;
            end
            
            if y<0
                y=0;
                first=1;
            elseif y>b
                y=b;
                first=1;
            end
        end
        
        
        
        boat_current = [x,y];
        
        distance = distance + dist(boat_last,boat_current);
        trajactory(t+1,:)=boat_current;
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
                    
                    for i5 = 1: numel(pts_in{i,i1})/2
                        boat_belief{i,i1}(pts_in{i,i1}(i5,2)+1,pts_in{i,i1}(i5,1)+1)=PDF(pts_in{i,i1}(i5,:),boat_current(i,:),range,sd); %prob_in{i3,i4}
                        boat_pass(pts_in{i,i1}(i5,2)+1,pts_in{i,i1}(i5,1)+1)=0;
                    end
                    belief{1,i1}=belief{1,i1}.*boat_belief{i,i1};
                    normalizer=1/sum(belief{1,i1}(:));
                    belief{1,i1}=belief{1,i1}.*normalizer;
                    if isnan(belief{1,i1})==1
                        belief{1,i1}=belief_last{1,i1};
                    end
                    belief_last{1,i1}=belief{1,i1};
                else
                    idx_in(i,i1)=0;
                    pts_in{i,i1}=overlapping(xy,boat_current(i,:),range+3*sd);
                    
                    for i5 = 1: numel(pts_in{i,i1})/2
                        boat_pass(pts_in{i,i1}(i5,2)+1,pts_in{i,i1}(i5,1)+1)=0;
                    end
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
            for i = 1:numel(c3)
                belief{1,c3(i)}=belief_last{1,c3(i)};
                variance(1,c3(i))=variance_last(1,c3(i));
            end
        elseif isempty(c4)==0
            for i=1:numel(c4)
                belief{1,c4(i)}=belief_last{1,c4(i)};
                variance(1,c4(i))=variance_last(1,c4(i));
            end
        end
        
        
        %%----------- verify change mode-------------------
        if enter_refine_mode ==0 && recal==0
            idx_delta=idx_in_last-idx_in;
            if numel(find(idx_delta>0))>0
                [r_r,c_r]=find((idx_delta)>0); 
                for i =1: numel(c_r)
                    if eta(c_r(i),t-2)<cl_level
                        enter_refine_mode=1;
                        first=1;
                    else
                        c_r(i)=[];
                    end
                end
            end
            idx_in_last = idx_in;    
        end
        
        if x==0 || x==b || y==0 || y==b && enter_refine_mode ==0
            recal =1;
        end

        variance_last = variance;
        
        %------------- estimate position ---------
        for  ii=1:s_number
            sumr=0;
            sumc=0;
            [r4,c4]=find(belief{1,ii}>belief_init);
            for i=1:numel(r4)
                sumr=sumr+belief{1,ii}(r4(i),c4(i))*(r4(i));
                sumc=sumc+belief{1,ii}(r4(i),c4(i))*(c4(i));
            end
            est_position(ii,2)=sumr-1;
            est_position(ii,1)=sumc-1;
        end
        
        %--------------- Confidence level ---------
        if isempty(r1)==0
            index=index+1;
            for ii1=1:s_number
                if numel(find(belief{1,ii1}>belief_init))==0
                    eta(ii1,t)=0;
                else
                    eta(ii1,t)=1-numel(find(belief{1,ii1}>belief_init))/ttgrid;
                end
                if t>1
                    if isnan(eta(ii1,t))==1
                        eta(ii1,t)=eta(ii1,t-1);
                    end
                end
                
            end
        else
            if index==1
                eta(:,t)=0;
            else
                eta(:,t)=eta(:,t-1);
            end
        end
        
        %%------------ detect parameter ---------
        dpmt=min(eta(:,t));
        if dpmt<0.3
            dpmt=0.3;
        elseif dpmt>0.85
            dpmt = 0.85;
        end 


        if isnan(position_r)==0
            if dist(boat_current,position_r)<10^-6
                rsearch_mode=1;
                start=1;
                first=1;
            end
        end
        
        %----------- Plot boat curretn position ------------------
        %     hFig = figure(1);
        %     set(gcf,'PaperPositionMode', 'auto')
        %     set(hFig,'Position',[100 200 1000 700])
        figure(1)
        subplot(2,2,1)
        colormap(jet)
        for i=1:b_number
            for ii=1:s_number
                if idx_in(i,ii)==1
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
        [anglemap,compass] = drawCompass(belief,boat_pass);
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
        drawPass(boat_pass);
        hold off
        zlim([0 1])
        axis tight
        title(['Probability Histogram after ',num2str(index-1),'th detection']);
        
        %     figure(4)
        subplot(2,2,4)
        for i2=1:s_number
            etapts = plot(t,eta(i2,t),'.','color',C(i2,:));
            hold on
            etaline = plot (1:t,eta(i2,1:t),'color',C(i2,:));
        end
        axis([0 500 0 1])
        title(['Confidence Level after ',num2str(index-1),'th detection'])
        %     delay between plots
        
        pause(.001)
        
        if numel(find(eta(:,t)>cl_level))==s_number  
            break
        end
    end
    auto=2;
end
while auto==2
    X=axis(joy, 1);     % X-axis is joystick axis 1
    Y=axis(joy, 2);     % Y-axis is joystick axis 2
    X_b = button(joy,1); % joystick X button to plot sensor position
    A_b = button(joy,2); % joystick A button to change speed
    B_b = button(joy,3); % joystick B button to plot trajactory
    Y_b = button(joy,4); % joystick Y button to change detect range and detect error
    if Y_b ==1
        if auto == 2
            auto =0;
        end
        %         range = input('Input detect range: ');
        %         sd = input('Input detect error: ');
        %         %-----calculate ttgrid --------
        %         tt=overlapping(xy,[b/2,b/2],range+3*sd);
        %         ttgrid = numel(tt)/2;
    end
    
    if A_b ==1
        v = input('Input speed: ');
    end
    
    %----------- Plot boat curretn position ------------------
    %     hFig = figure(1);
    %     set(gcf,'PaperPositionMode', 'auto')
    %     set(hFig,'Position',[100 200 1000 700])
    figure(1)
    subplot(2,2,1)
    colormap(jet)
    %         for i=1:b_number
    %             for ii=1:s_number
    %                 if idx_in(i,ii)==1
    %                     plot(pts_in{i,ii}(:,1),pts_in{i,ii}(:,2),'r.');
    %                     hold on
    %                 end
    %             end
    %         end
    plot(boat_current(1,1),boat_current(1,2),'>','MarkerSize',10,'MarkerFaceColor','y');
    hold on
    if X_b == 1
        plot(sensor_position(:,1),sensor_position(:,2),'b+');
    end
    if B_b ==1
        plot(trajactory(:,1),trajactory(:,2),'b-');
    end
    [anglemap,compass] = drawCompass(belief,boat_pass);
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
    drawPass(boat_pass);
    hold off
    zlim([0 1])
    axis tight
    title(['Probability Histogram after ',num2str(index-1),'th detection']);
    
    %     figure(4)
    subplot(2,2,4)
    for i2=1:s_number
        etapts = plot(t,eta(i2,t),'.','color',C(i2,:));
        hold on
        etaline = plot (1:t,eta(i2,1:t),'color',C(i2,:));
    end
    axis([0 2000/v 0 1])
    title(['Confidence Level after ',num2str(index-1),'th detection'])
    %     delay between plots
    
    pause(.001)
end