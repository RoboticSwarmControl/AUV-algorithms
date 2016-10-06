%% distance calculation
function d = dist(a,b)
       d = sqrt((a(1,1)-b(1,1))^2+(a(1,2)-b(1,2))^2);
end