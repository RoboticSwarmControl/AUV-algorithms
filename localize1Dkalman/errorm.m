function dist = errorm(sensor,boat)
l = sqrt((sensor(1,1)-boat(1,1))^2);
dist = l+randn(1,1)*5;
end