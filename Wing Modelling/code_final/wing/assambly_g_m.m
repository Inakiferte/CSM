function [K,M,l,R,Me,w,Nm,Ba_p,Bs_p,Bt_p,Bb_p,Ka,Kb,Ks,Kt] = assambly_g_m(N_dof,n_elem,xn,Tn,y_p,E,A,G,rho,J,Iz,Iy,ky,kz,kt)
% Assembly the global matrices
K     = sparse(N_dof,N_dof);
M     = sparse(N_dof,N_dof);
for e=1: n_elem
    % a) Compute rotation matrix
        l(e) = norm(xn(Tn(e,2),:) - xn(Tn(e,1),:));
        i_p = (transpose(xn(Tn(e,2),:)) - transpose(xn(Tn(e,1),:))) / l(e);
        j_p = y_p;
        k_p = cross(i_p,j_p);

        R_p = zeros(6,6);
        R_p(1,1:3) = [i_p];
        R_p(2,1:3) = [j_p];
        R_p(3,1:3) = [k_p];
        R_p(4,4:6) = [i_p];
        R_p(5,4:6) = [j_p];
        R_p(6,4:6) = [k_p];
        R_p = transpose(R_p);

        R(1:6,1:6,e) = R_p;
        R(7:12,7:12,e) = R_p;
    % b) Compute shape function derivatives
        N_x(1) = -1 / l(e);
        N_x(2) = 1 / l(e);
    % C) Compute each element matrix
        % c1) Axial component of stiffness matrix
        Ba_p(1,:,e) = [N_x(1),0,0,0,0,0,N_x(2),0,0,0,0,0];
        Ca_p = E * A;
        Ka(:,:,e) = l(e) * transpose(R(:,:,e)) * transpose(Ba_p(1,:,e)) * Ca_p * Ba_p(1,:,e) * R(:,:,e);
        % c2) Bending component of stiffness matrix
        Bb_p(:,:,e) = [0,0,0,0,N_x(1),0,0,0,0,0,N_x(2),0;
                       0,0,0,0,0,N_x(1),0,0,0,0,0,N_x(2)];
        Cb_p(1,1) = Iy;
        Cb_p(1,2) = 0;
        Cb_p(2,1) = 0;
        Cb_p(2,2) = Iz;
        Cb_p = E * Cb_p;
        Kb(:,:,e) = l(e) * transpose(R(:,:,e)) * transpose(Bb_p(:,:,e)) * Cb_p * Bb_p(:,:,e) * R(:,:,e);
        % c3) Shear component of stiffness matrix
        N = 0.5;
        Bs_p(:,:,e) = [0,N_x(1),0,0,0,-N,0,N_x(2),0,0,0,-N;
                       0,0,N_x(1),0,N,0,0,0,N_x(2),0,N,0];
        Cs_p(1,1) = ky;
        Cs_p(1,2) = 0;
        Cs_p(2,1) = 0;
        Cs_p(2,2) = kz;
        Cs_p = G * A * Cs_p;
        Ks(:,:,e) = l(e) * transpose(R(:,:,e)) * transpose(Bs_p(:,:,e)) * Cs_p * Bs_p(:,:,e) * R(:,:,e);
        % c4) Torsion component of stiffness matrix
        Bt_p(1,:,e) = [0,0,0,N_x(1),0,0,0,0,0,N_x(2),0,0]; 
        Ct_p = G * J * kt;
        Kt(:,:,e) = l(e) * transpose(R(:,:,e)) * transpose(Bt_p(1,:,e)) * Ct_p * Bt_p(1,:,e) * R(:,:,e);
        % c5) Mass matrix
        xhi = [-1 / sqrt(3);1 / sqrt(3)];
        w   = [1;1];
        rho_p = zeros(6,6);
        rho_p(1,1:6) = [A,0,0,0,0,0];
        rho_p(2,1:6) = [0,A,0,0,0,0];
        rho_p(3,1:6) = [0,0,A,0,0,0];
        rho_p(4,1:6) = [0,0,0,J,0,0];
        rho_p(5,1:6) = [0,0,0,0,Iy,0];
        rho_p(6,1:6) = [0,0,0,0,0,Iz];
        rho_p = rho * rho_p;
        Me(:,:,e) = zeros(12,12);
        for k=1:2
            Ng(1) = (1 - xhi(k)) / 2;
            Ng(2) = (1 + xhi(k)) / 2;
            Ident = eye(6,6);
            Nm(:,:,e,k) = [Ng(1)*Ident,Ng(2)*Ident];
            Me(:,:,e) = Me(:,:,e) + w(k) * l(e) * transpose(R(:,:,e)) * transpose(Nm(:,:,e,k)) * rho_p * Nm(:,:,e,k) * R(:,:,e) / 2;
        end
        % d) Assembly to global matrices
        for j=1:6
            Idof(j,1) = 6*(Tn(e,1) - 1) + j;
            Idof(6 + j,1) = 6*(Tn(e,2) - 1) + j;
        end
        K(Idof(:,1),Idof(:,1)) = K(Idof(:,1),Idof(:,1)) + Ka(:,:,e) + Kb(:,:,e) + Ks(:,:,e) + Kt(:,:,e);
        M(Idof,Idof) = M(Idof,Idof) + Me(:,:,e);
end