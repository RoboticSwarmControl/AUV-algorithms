function hHistogram = updthist(s_number,belief,hHistogram)
for i = 1:s_number
    if isempty(find(belief{1,i})>0)==0
        set(hHistogram{1,i},'xdata',belief{1,i});
    end
    hold on
end
end