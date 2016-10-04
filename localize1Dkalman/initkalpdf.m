function [hKalpdf,hKalpts]=initkalpdf(n,mean_new,cov_new)
hKalpdf = zeros(n,1);
hKalpts = zeros(n,1);
for i=1:n
    range= mean_new(i,1)-40:.1:mean_new(i,1)+40;
    pdf = exp(- 0.5 * ((range - mean_new(i,1)) / cov(i,1)) .^ 2) / (cov_new(i,1) * sqrt(2 * pi));
    hKalpdf(i) = plot(range,pdf,'k');
    hKalpts(i) = plot(mean(i,1),0,'k+');
end