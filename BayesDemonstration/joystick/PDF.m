function prob = PDF(pts_in,boat_current,range,sd)
d=dist(pts_in,boat_current);
if d<=30
    prob = exp(-(range-range)^2/(2*sd^2))/(sqrt(2*sd^2*pi));
else
    prob= exp(-(d-range)^2/(2*sd^2))/(sqrt(2*sd^2*pi));
end
end