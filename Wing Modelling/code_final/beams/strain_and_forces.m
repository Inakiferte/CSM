function [epsilon_a,epsilon_s,epsilon_t,epsilon_b,Fx,Fy,Fz,Mx,My,Mz] = strain_and_forces(n_elem,Tn,u_hat,R,Ba_p,Bs_p,Bt_p,Bb_p,Ka,Kb,Ks,Kt)
for e=1: n_elem
    % a) Get element displacements
    for j=1:6
        Idof(j,1) = 6*(Tn(e,1) - 1) + j;
        Idof(6 + j,1) = 6*(Tn(e,2) - 1) + j;
    end
    u_hate = u_hat(Idof,1);
    % b) Get each strain component
    epsilon_a(1,e) = Ba_p(1,:,e) * R(:,:,e) * u_hate;
    epsilon_s(:,e) = Bs_p(:,:,e) * R(:,:,e) * u_hate;
    epsilon_t(1,e) = Bt_p(1,:,e) * R(:,:,e) * u_hate;
    epsilon_b(:,e) = Bb_p(:,:,e) * R(:,:,e) * u_hate;
    % c) Get internal forces and moments at each element node
    f_hatint(:,e) = R(:,:,e) * (Ka(:,:,e) + Kb(:,:,e) + Ks(:,:,e) + Kt(:,:,e)) * u_hate;
    Fx(:,e) = [f_hatint(1,e); -f_hatint(7,e)];
    Fy(:,e) = [f_hatint(2,e); -f_hatint(8,e)];
    Fz(:,e) = [f_hatint(3,e); -f_hatint(9,e)];
    Mx(:,e) = [f_hatint(4,e); -f_hatint(10,e)];
    My(:,e) = [f_hatint(5,e); -f_hatint(11,e)];
    Mz(:,e) = [f_hatint(6,e); -f_hatint(12,e)];
end
end