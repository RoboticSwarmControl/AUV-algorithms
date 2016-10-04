function [hEstpdf,hEstpts]=updmpdf(n,mean,cov,hEstpdf,hEstpts)
for i=1:n
    range= mean(i)-40:.1:mean(i)+40;
    pdf = exp(- 0.5 * ((range - mean(i)) / cov(i)) .^ 2) / (cov(i) * sqrt(2 * pi));
    set (hEstpdf(i),'xdata',range,'ydata',pdf);
    set (hEstpts(i),'xdata',mean(i),'ydata',0);
end
end