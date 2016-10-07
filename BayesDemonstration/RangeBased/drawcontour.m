function drawcontour(belief)
contour(belief,5);
xtick = 0:10:100;
xticklabels = {'-50','-40','-30','-20','-10','0','10','20','30','40','50'};
set(gca,'XTick', xtick);
set(gca,'XTickLabel', xticklabels);
ytick =  0:10:100;
yticklabels = {'-50','-40','-30','-20','-10','0','10','20','30','40','50'};
set(gca,'YTick', ytick);
set(gca,'YTickLabel', yticklabels);
axis equal
end