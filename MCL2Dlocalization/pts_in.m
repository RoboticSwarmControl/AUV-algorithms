%% point in overlapping area
    function xy= pts_in(dist,cent)
        [~,c]= find(dist==min(dist(:)));
        [~,col]=size(dist);
        d=min(dist);
        range_min=cent(c,:)-d;
        if range_min(1) < 0
            range_min(1) = 0;
        end
        if range_min(2) < 0
            range_min(2) = 0;
        end
        range_max=cent(c,:)+d;
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
       xy=overlapping(xy,cent(c,:),d);
       dist_new = dist;
       dist_new(:,c)=[];
       cent_new(:,:) = cent(:,:);
       cent_new(c,:)=[];
       for i=1:col-1
           r=dist_new(1,i);
       xy=overlapping(xy,cent_new(i,:),r);
       end
       
    end