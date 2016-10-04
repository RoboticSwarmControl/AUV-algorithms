function hWorldFig= drawWorld(belief)
figure(2)
colormap(jet)
hWorldFig =bar3(belief,1);
set(gca,'YDir','normal');
[obsy,obsx] = find(belief == 0);
for i = 1:numel(obsy)
   rectangle('Position',[obsx(i)-0.5,obsy(i)+0.5,1,1],'FaceColor','b')
end
end
