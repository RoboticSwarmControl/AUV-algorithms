function hSppts = updSppts(s_number,xy,hSppts)
       
    for i=1:s_number
       set(hSppts(i), 'xdata', xy(:,1),'ydata',xy(:,2));
    end
end