function [hDfpts,hDfline]=upddfdist(n,t,diff,hDfpts,hDfline)
if t==0
    col =1;
else
    col = t/0.1+1;
end
col=round(col);
for i=1:n
    set (hDfpts(i),'xdata',t,'ydata',diff(i,col));
    set(hDfline(i),'xdata',0:.1:t,'ydata',diff(i,1:col));
end
 
end