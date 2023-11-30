function [output] = body_force(x,c,lambda)
output = c * exp(lambda * x);
end