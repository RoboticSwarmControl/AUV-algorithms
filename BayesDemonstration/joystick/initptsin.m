function hPtsin = initptsin(s_number,b_number,pts_in)
hPtsin = zeros(b_number,s_number);
for i=1:b_number
    for ii=1:s_number
        if isempty(find(pts_in{i,ii})>0)==0
            hPtsin(i,ii)=plot(pts_in{i,ii}(:,1),pts_in{i,ii}(:,2),'r.');
            hold on
        end
    end
end
end