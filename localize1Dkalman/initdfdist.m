function [hDfpts,hDfline]=initdfdist(n,t,diff)
if t==0
    col =1;
else
    col = t/0.1+1;
end
col=round(col);
hDfpts = zeros(n,1);
hDfline = zeros(n,1);
for i = 1:n
    hDfpts(i) = plot (t,diff(i,col),'.');
    hDfline(i) = plot (0:.1:t,diff(i,1:col));
end

end

