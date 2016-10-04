function [mu_new,cov_new]=kalmanfilter(mu,cov_e,mum)
v=mum-mu;
s=cov_e+30;
r=cov_e/s;
mu_new = mu+r*v;
cov_new = cov_e -r*cov_e;
end
