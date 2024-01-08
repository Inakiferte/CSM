function [U_s] = modal_order_reduction(Im,Ndof,f,lambda,phi)
% 9.2 Project the system onto the selected set of nodes
Nw = 1;
U_s = zeros(Ndof,Nw);
for k=1: size(f,2)
    for j=1:length(Im)
        alpha(j,k) = phi(:,Im(j))' * f(:,k) / lambda(j);
        U_s(:,k)   = U_s(:,k) + (phi(:,Im(j)) * alpha(j,k));
    end
end
end