function [hKalpdf,hKalpts]=updkalpdf(n,mean_new,cov_new,hKalpdf,hKalpts)
for i=1:n
    range= mean_new(i)-40:.1:mean_new(i)+40;
    pdf = exp(- 0.5 * ((range - mean_new(i)) / cov_new(i)) .^ 2) / (cov_new(i) * sqrt(2 * pi));
    set (hKalpdf(i),'xdata',range,'ydata',pdf);
    set (hKalpts(i),'xdata',mean_new(i),'ydata',0);
end
end