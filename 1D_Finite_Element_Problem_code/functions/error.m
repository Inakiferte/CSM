function [output1, output2] = error(u, u_ana, s, s_ana)
output1 = abs(u(end)-u_ana(end))/abs(u(end));
output2 = abs(s(end)-s_ana(end))/abs(s(end));
end