function hCline = updtcontour(s_number,belief,hCline)
for ii=1:s_number
    set(hCline{ii,1},'xdata',belief{1,ii},'ydata',5);
end
end