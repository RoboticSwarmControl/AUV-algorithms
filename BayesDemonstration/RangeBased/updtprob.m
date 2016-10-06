function prob = updtprob(d_c,d_m,sd)
prob = exp(-(d_c-d_m)^2/(2*sd^2))/sqrt(2*sd^2*pi);
end