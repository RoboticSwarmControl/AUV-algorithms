function [confPts, confLine] = updtconf(s_number,index,eta,confPts,confLine)
for i2=1:s_number
    set(confPts(i2),'xdata',index,'ydata',eta(i2,index));
    hold on
    set(confLine(i2),'xdata',1:index,'ydata',eta(i2,1:index));
end
end