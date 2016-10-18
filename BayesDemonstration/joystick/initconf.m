function [confPts, confLine] = initconf(s_number,index,eta)
C = hsv(s_number);
confPts = zeros(s_number,1);
confLine = zeros(s_number,1);
    for i2=1:s_number
        confPts(i2) = plot(index,eta(i2,index),'.','color',C(i2,:));
        hold on 
        confLine(i2) = plot (1:index,eta(i2,1:index),'color',C(i2,:));
    end
end