function [hMspdf,hMspts]=initmpdf(n,mean,cov)
hMspdf = zeros(n,1);
hMspts = zeros(n,1);
for i=1:n
    range= mean(i,1)-40:.1:mean(i,1)+40;
    pdf = exp(- 0.5 * ((range - mean(i,1)) / cov(i,1)) .^ 2) / (cov(i,1) * sqrt(2 * pi));
    hMspdf(i) = plot(range,pdf,'r');
    hMspts(i) = plot(mean(i,1),0,'k+');
end