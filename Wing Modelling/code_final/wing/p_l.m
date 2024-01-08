function [p] = p_l(x,y)
% Input values
alpha = deg2rad(10);
c     = 0.8;
p_x = px(x);

% Output values
p = 4 * alpha * p_x * ((1 - (y / c))^4 - (1 / 4) * sqrt(1 - (y / c)));

end