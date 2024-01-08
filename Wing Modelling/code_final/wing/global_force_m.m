function [f_hat] = global_force_m(N_dof,Fe,n,Be,n_elem,Qe,Me,w,l,R,Nm,Tn)
% 3.1 Points loads
f_hat = zeros(N_dof,1);
for q=1:size(Fe)
    f_hat(6*(Fe(q,2) - 1) + Fe(q,3),1) = f_hat(6*(Fe(q,2) - 1) + Fe(q,3),1) + Fe(q,1);
end
% 3.2 Nodal distributed forces
Q = zeros(n,6);
% No distributed load matrix
for r=1:size(Qe,1)
     Q(Qe(r,2),Qe(r,3)) = Q(Qe(r,2),Qe(r,3)) + Qe(r,1);
end

% 3.3 Nodal body forces
B = zeros(n,6);
for s=1:size(Be,1)
    B(Be(s,2),Be(s,3)) = B(Be(s,2),Be(s,3)) + Be(s,1);
end

% 3.4 Assembly process
for e=1:n_elem
    % a) Compute element force vector
    b(:,e) = transpose([B(Tn(e,1),:), B(Tn(e,2),:)]);
    qp(:,e) = transpose([Q(Tn(e,1),:), Q(Tn(e,2),:)]); 
    fe_hat(:,e) = Me(:,:,e) * b(:,e);
    for k=1:2
        fe_hat(:,e) = fe_hat(:,e) + w(k) * l(e) * transpose(R(:,:,e)) * transpose(Nm(:,:,e,k)) * Nm(:,:,e,k) * R(:,:,e) * qp(:,e) / 2;
    end
    % b) Assembly to globar force vector
    for j=1:6
        Idof(j,1) = 6*(Tn(e,1) - 1) + j;
        Idof(6 + j,1) = 6*(Tn(e,2) - 1) + j;
    end
    f_hat(Idof,1) = f_hat(Idof,1) + fe_hat(:,e);
end

end