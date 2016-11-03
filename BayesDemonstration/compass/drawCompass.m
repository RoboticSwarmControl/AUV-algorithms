function [anglemap1,compass] = drawCompass(belief,boat_pass)
n=numel(belief);
[r,c]=size(belief{1,1});
r_c=(r-1)/2+1;
c_c=(c-1)/2+1;
r1=cell(1,72);
c1=cell(1,72);
anglemap=zeros(size(r,c));
% prob = cell(1,n);
% pass_idx=zeros(size(boat_pass));
% prob_pass=zeros(2,72);
Z = zeros(2,72);

    for i1 = 1:r
        for i2=1:c
            s=atan2(i2-1-c_c,i1-1-r_c)/pi*180;
            if s<0
                s=s+360;
            end
            anglemap(i1,i2)=ceil(s/5);
        end
    end
    anglemap1=anglemap;
    
    for i = 1:72
        [r1{1,i},c1{1,i}]=find(anglemap==i);    
    end
    
    for i=1:n
        for i1=1:72
            num=numel(r1{1,i1});
            for i2=1:num
            Z(1,i1)=Z(1,i1)+belief{1,i}(r1{1,i1}(i2),c1{1,i1}(i2));
            Z(2,i1)=Z(1,i1);
            end
        end
    end
    
    for i =1:72
        num=numel(r1{1,i});
        for i1 =1:num
            r2=r1{1,i}(i1);
            c2=c1{1,i}(i1);
            Z(1,i)=Z(1,i)+boat_pass(c2,r2);
            Z(2,i)=Z(1,i);
        end
    end
    
    

th = linspace(0,2*pi,72);
r = r_c-7:5:r_c-2;
[TH,R] = meshgrid(th,r);
[X,Y] = pol2cart(TH,R);
% xrange = 0:200;
% yrange = 0:200;
% [X1,Y1] = ndgrid(xrange,yrange);
% Z1=boat_pass;
% background = Z1;
% for i= 1:n
%     Z1=Z1+belief{1,i};
% end


%  ax1 = axes;
% bPlot=surf(Y1,X1,Z1);
% view(2)
% set(bPlot,'FaceColor','interp','EdgeColor','interp');
% colormap(ax1,'pink');
% % cb1 = colorbar(ax1,'Location','westoutside');
% % for ks = 1:length(bPlot)
% %     zdata = get(bPlot(ks),'ZData');
% %     set(bPlot(ks),'CData', zdata);
% %     set(bPlot(ks),'FaceColor', 'interp');
% % end

% ax2 = axes;
compass = Z;
Z=Z./sum(Z(1,:));
colormap(jet)
cPlot=surf(Y+r_c-1,X+r_c-1,Z);
view(2)
set(cPlot,'FaceColor','interp','EdgeColor','interp');
% colormap(ax2,'jet');
%  cb2 = colorbar(ax2,'Location','eastoutside');
for ks = 1:length(cPlot)
    zdata = get(cPlot(ks),'ZData');
    set(cPlot(ks),'CData', zdata);
    set(cPlot(ks),'FaceColor', 'interp');
end

% linkaxes([ax1,ax2])

% ax2.Visible = 'off';
% ax2.XTick = [];
% ax2.YTick = [];


axis square
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);
% axis([ax1,ax2],'square')
end



