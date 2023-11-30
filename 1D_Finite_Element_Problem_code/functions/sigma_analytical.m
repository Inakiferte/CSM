function [sigma] = sigma_analytical(x,c,A,lambda,L)
    sigma = (-c / A) * ((1 / lambda) * exp(lambda * x) - (exp(lambda * L) / lambda));
end