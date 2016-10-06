%function BayesDemonstrationRangefree()
%BayesDemonstration is used to demonstrate Bayes' rule for an underwater
%localization sestup
%
%The algorithm uses range-free localization
%
%
% Subfunctions: 
%
% Author: Haoran Zhao, Master, Electrical Engineering
% University of Houston
% email address: zhaohaorandl@gmail.com
% Last revision: 4-Oct-2016

%------------- BEGIN CODE --------------
clear
clc
format compact
%------------- Main variable -----------
s_number = 2; %sensor number
b_number = 3; %boat number

%generate sensor position
a=0;
b=100;
c=20;
d=80;
x1=c+(d-c)*rand(s_number,1);
y1=c+(d-c)*rand(s_number,1);
sensor_position=[x1(:),y1(:)];

%generate map area
xrange = 0:100;
yrange = 0:100;
[X,Y]=ndgrid(xrange,yrange);
xy=[X(:),Y(:)];

range = 30; %boat sensor detect range
sd = 1; %sensor error 

belief = cell(1,s_number); % Bayes' rule belief
belief_init = ones(101,101)./numel(ones(101,101));

for i=1:s_number
    belief{1,i}=belief_init;
end

boat_belief = cell(b_number,s_number); % Boat knows if there is sensor in the range
b_position = zeros(b_number,s_number); % Boat history position 
pts_in = cell(b_number,s_number); % find the possible position in the range
prob_in = cell (b_number,s_number);
idx_in = zeros(b_number,s_number); % record if sensor in range with 0 or 1

distance_bws=zeros(b_number,s_number);
variance = zeros(1,s_number);
est_postion = zeros(s_number,2);
%------------ Main Code ------------
for t=1:1000
    %generate boat position
    x2=a+(b-a)*rand(b_number,1);
    y2=a+(b-a)*rand(b_number,1);
    boat_position = [x2(:),y2(:)];
    
    pts_in=cell(b_number,s_number);
    belief_last=cell(1,s_number);
    
    for i = 1:b_number
        for i1 = 1:s_number
            %calculate distance bwt boat and sensor
            distance_bws(i,i1) = dist(sensor_position(i1,:),boat_position(i,:));
            detectrange = range+randn(1,1)*sd;
            if distance_bws(i,i1)< detectrange || distance_bws(i,i1) == detectrange
                idx_in(i,i1)=1;
            else
                idx_in(i,i1)=0;
            end
        end
    end
    
    [r1,c1]=find(idx_in==1);
    if isempty(r1)~=1
        for i2=1:numel(r1)
            %find possible position in the range
            pts_in{r1(i2),c1(i2)}=overlapping(xy,boat_position(r1(i2),:),range);  
        end
        for i3=1:b_number
            for i4=1:s_number
                if isempty(pts_in{i3,i4})==0
                 prob_in{i3,i4}=2/numel(pts_in{i3,i4});
                 boat_belief{i3,i4}=zeros(101,101);
                for i5 = 1:numel(pts_in{i3,i4})/2
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
    
    for i6=1:s_number
    [r2,c2]=find(belief{1,i6}>0);
    variance(1,i6)=var(r2(:),c2(:));
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
    
%---------------plot image ---------------
figure(1)
plot (boat_position(:,1),boat_position(:,2),'>','markerfacecolor','y','markersize',10);
hold on 
plot (sensor_position(:,1),sensor_position(:,2),'b+');
% for i8=1:b_number
%     for i9=1:s_number
%         if isempty(find(pts_in{i8,i9})>0)==0
%             plot(pts_in{i8,i9}(:,1),pts_in{i8,i9}(:,2),'r.');
%         end
%     end
% end
for ii=1:s_number
    contour(belief{1,ii},5);
end
hold off
axis equal
axis([0 100 0 100])
legend('Boat position','sensor position','Estimated Area');
title({['Histogram after ',num2str(t),'  range detection'],['Histogram variance = ', num2str(variance)],})


figure(2)
for i7=1:s_number
    if isempty(find(belief{1,i7})>0)==0
    redrawWorld(belief{1,i7});
    end
    hold on
end
hold off
zlim([0 1])
axis tight
title({['Histogram after ',num2str(t),'  range detection'],['Histogram variance = ', num2str(variance)],})

figure(3)
for i7=1:s_number
    if isempty(find(belief{1,i7})>0)==0
    redrawWorlds(belief{1,i7});
    end
    hold on
end
hold off
zlim([0 1])
axis tight
title({['Histogram after ',num2str(t),'  range detection'],['Histogram variance = ', num2str(variance)],})

for ii=1:s_number
    disp(['Sensor position : (x,y)= [',num2str(sensor_position(ii,1)),',',num2str(sensor_position(ii,2)),']']) 
    disp(['Estimate position after ',num2str(t),' range detection: (x,y)= [',num2str(est_position(ii,1)),',',num2str(est_position(ii,2)),']'])
end
pause(2);
end


    
    
            








