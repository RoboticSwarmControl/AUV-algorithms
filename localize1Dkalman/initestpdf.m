function [hEstpdf,hEstpts]=initestpdf(n,mean,cov)
hEstpdf = zeros(n,1);
hEstpts = zeros(n,1);
for i=1:n
    range= mean(i,1)-100:.1:mean(i,1)+100;
    pdf = exp(- 0.5 * ((range - mean(i,1)) / cov(i,1)) .^ 2) / (cov(i,1) * sqrt(2 * pi));
    hEstpdf(i) = plot(range,pdf,'b');
    hEstpts(i) = plot(mean(i,1),0,'b+');
end
