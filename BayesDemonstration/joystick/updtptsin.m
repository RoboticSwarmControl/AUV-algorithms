function hPtsin = updtptsin(s_number,b_number,pts_in,hPtsin)
for i=1:b_number
    for ii=1:s_number
        if isempty(find(pts_in{i,ii})>0)==0
            set(hPtsin(i,ii),'xdata',pts_in{i,ii}(:,1),'ydata',pts_in{i,ii}(:,2));
            hold on
        end
    end
end
end