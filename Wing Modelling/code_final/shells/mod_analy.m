function [omega,phi,lambda] = mod_analy(Nm,K,M,If,N_dof)
% Modal analysis function

% Tip: Compute simmetric K and M
K = (K + transpose(K)) / 2;
M = (M + transpose(M)) / 2;

% 7.1 Solve the eingenvalue and eigenvector problem, diagonalize
[V,D] = eigs(K(If,If),M(If,If),Nm,'sm');

% 7.2 Obtain squared natural frequencies
phi    = zeros(N_dof,Nm);
lambda = zeros(1,Nm);

for k=1:size(V,2)
    % Mass-normalized vibration mode
    phi(If,k) = V(:,k) / sqrt(transpose(V(:,k)) * M(If,If) * V(:,k));
    % Squared natural frequency
    lambda(k) = D(k,k);
end

% Obtain mode frecuencies
omega = sqrt(lambda);

end