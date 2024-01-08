function [p] = px(x)
p_infty = 1.25 * 1e4;             % [Pa]
b = 5.5;                          % NACA0012 value

if x<=(0.8 * b)
    p = p_infty;
else
    p = p_infty * cos((pi / 2) * ((x - 0.8 * b) / (0.2 * b) ));
end
end