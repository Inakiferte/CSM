function [u] = u_analytical(x,c,E,A,lambda,L)
    u = (-c / (E * A)) * ((1 / lambda^2) * exp(lambda * x) - (x / lambda) * exp(lambda * L) - (1 / lambda^2));
end