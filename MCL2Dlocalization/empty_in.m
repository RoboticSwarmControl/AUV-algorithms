function xy = empty_in(xy,distance_error,A_position)
if isempty(xy)==1
       [~,col]=size(distance_error);
       [~,c]= find(distance_error==min(distance_error(:)));
       dist_new = distance_error;
       dist_new(:,c)=[];
       cent_new = A_position;
       cent_new(c,:)=[];
       [~,c_new]= find(dist_new==min(dist_new(:)));
        d=min(dist_new);
        range_min=cent_new(c_new,:)-d;
        if range_min(1) < 0
            range_min(1) = 0;
        end
        if range_min(2) < 0
            range_min(2) = 0;
        end
        range_max=cent_new(c_new,:)+d;
        if range_max(1) > 100
            range_max(1) = 100;
        end
        if range_max(2) > 100
            range_max(2) = 100;
        end
        xrange = ceil(range_min(1)):floor(range_max(1));
        yrange = ceil(range_min(2)):floor(range_max(2));
        [X,Y]=ndgrid(xrange,yrange);
        xy = [X(:),Y(:)];
        xy=overlapping(xy,A_position(c_new,:),d);
       dist_new(:,c_new)=[];
       cent_new(c_new,:)=[];
       for j=1:col-2
           r=dist_new(1,j);
       xy=overlapping(xy,cent_new(j,:),r);
       end       
end

end