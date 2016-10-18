function hHistogram = inithist(s_number,belief)
hHistogram =cell(1,s_number);
for i = 1:s_number
    if isempty(find(belief{1,i})>0)==0
        hHistogram{1,i}= redrawWorlds(belief{1,i});
    end
    hold on
end
end