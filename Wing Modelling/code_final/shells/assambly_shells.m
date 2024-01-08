function [K,M,R,Me,NN,Bs_p,Bb_p,Kb,Ks,S_4,Bmn_p,Bmt_p] = assambly_shells(Ndof,n_elem,xn,Tn,Tm,E,rho,h1,h2,h3,nu)
% Obtain system matrices
K = sparse(Ndof,Ndof);
M = sparse(Ndof,Ndof);
ze = zeros(3,2);
ze_p = zeros(3,3);
ze_pp = zeros(5,6);
for e = 1:n_elem
    % Assembly process:
    % a) Compute Rotational Matrix
    S = cross((xn(Tn(e,3),:)' - xn(Tn(e,1),:)'),(xn(Tn(e,4),:)'-xn(Tn(e,2),:)'))/2;
    k_p = S/norm(S);                    % Normal vector of flat shell element
    d = (xn(Tn(e,2),:)' + xn(Tn(e,3),:)' - xn(Tn(e,4),:)' - xn(Tn(e,1),:)')/2;
    i_p = d/norm(d);
    j_p = cross(k_p,i_p);
    R_p = [i_p j_p k_p ze ; ze_p i_p j_p]';
    R(:,:,e) = [R_p ze_pp ze_pp ze_pp ; ze_pp R_p ze_pp ze_pp ; ze_pp ze_pp R_p ze_pp ; ze_pp ze_pp ze_pp R_p];

    % b) get nodal coefficients for the shape functions
    a = [-1 1 1 -1];
    b = [-1 -1 1 1];

    % c) Compute element matrices
    % c1) 1. Gauss point quadrature matrices
    N1 = [1, 1, 1, 1]'/4;
    N1_xhi = a/4;
    N1_eta = b/4;
    J_1 = zeros (2,2);

    for i=1:4
        J_1 = J_1 + [N1_xhi(i);N1_eta(i)]*xn(Tn(e,i),:)*[i_p j_p];
    end
    N1_xp = inv(J_1)*[N1_xhi;N1_eta];
    S1 = 4*det(J_1);           % Area associated to Gauss point

    % c1.1 Shear component of stiffness matrix
    for i=1:4
        bs_p(:,:,i) = [0 0 N1_xp(1,i) 0 N1(i); 0 0 N1_xp(2,i) -N1(i) 0];
    end
    if Tm(e) ==1
        h = h1;
    elseif Tm(e) == 2
        h = h2;
    else 
        h = h3;
    end
    Cs_p = [1 0; 0 1] * 5 * h * E/(12*(1 + nu));
    Bs_p(:,:,e) = [bs_p(:,:,1) bs_p(:,:,2) bs_p(:,:,3) bs_p(:,:,4)];
    Ks(:,:,e) = S1 * R(:,:,e)' * Bs_p(:,:,e)' * Cs_p * Bs_p(:,:,e) * R(:,:,e);

    % c1.2) Mmebrane transverse component of stiffness matrix
    for i=1:4
        bmt_p(:,:,i) = [N1_xp(2,i) N1_xp(1,i) 0 0 0];
    end
    Cmt_p = h * E/(2*(1+nu));
    Bmt_p(:,:,e) = [bmt_p(:,:,1) bmt_p(:,:,2) bmt_p(:,:,3) bmt_p(:,:,4)];
    Km(:,:,e) = S1 * R(:,:,e)' * Bmt_p(:,:,e)' * Cmt_p * Bmt_p(:,:,e) * R(:,:,e);

    % c2) 4 Gauss points quadrature matrices
    Kb (:,:,e) = zeros(24,24);
    Me (:,:,e) = zeros(24,24);
    xhi_4 = [-1 1 1 -1]/sqrt(3);
    eta_4 = [-1 -1 1 1]/sqrt(3);
    w_4 = [1 1 1 1];
    for k=1:4
        J_4 = zeros (2,2);
        for i=1:4
            N_4(i) = (1+a(i)*xhi_4(k))*(1+b(i)*eta_4(k))/4;
            N_4_xhi(1,i) = a(i)*(1+b(i)*eta_4(k))/4;
            N_4_eta(1,i) = b(i)*(1 + a(i)*xhi_4(k))/4;
            J_4 = J_4 + [N_4_xhi(i); N_4_eta(i)]*(xn(Tn(e,i),:))* [i_p j_p];
        end
        N_4_xp = inv(J_4)*[N_4_xhi;N_4_eta];
        S_4(e,k) = w_4(k)*det(J_4);        % Area associated to Gauss point

        % c2.1) Membrane normal component of stiffness matrix

        for i=1:4
            bmn_p(:,:,i) = [N_4_xp(1,i) 0 0 0 0; 0 N_4_xp(2,i) 0 0 0];
        end
        Cmn_p = [1 nu; nu 1] * h * E/(1 - nu^2);
        Bmn_p (:,:,e,k) = [bmn_p(:,:,1) bmn_p(:,:,2) bmn_p(:,:,3) bmn_p(:,:,4)];
        Km (:,:,e) = Km(:,:,e) + S_4(e,k) * R(:,:,e)' * Bmn_p(:,:,e,k)' * Cmn_p * Bmn_p(:,:,e,k) * R(:,:,e);

        % c2.2) Bending component of stiffness matrix
        for i=1:4
            bb_p(:,:,i) = [0 0 0 0 N_4_xp(1,i); 0 0 0 N_4_xp(2,i) 0; 0 0 0 -N_4_xp(1,i) N_4_xp(2,i)];
        end
        
        Cb_p = [1 nu 0; nu 1 0; 0 0 (1-nu)/2] * h^3 * E/ (12 *(1-nu^2));
        Bb_p(:,:,e,k) = [bb_p(:,:,1) bb_p(:,:,2) bb_p(:,:,3) bb_p(:,:,4)];
        Kb(:,:,e) = Kb(:,:,e) + S_4(e,k) * R(:,:,e)' * Bb_p(:,:,e,k)' * Cb_p * Bb_p(:,:,e,k) * R(:,:,e);

        % c2.3) Mass matrix

        for i=1:4
            N(:,:,i) = N_4(i) * eye(5); 
        end
        rho_p = rho * h * [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 (h^2)/12 0; 0 0 0 0 (h^2)/12];
        NN(:,:,e,k) = [N(:,:,1) N(:,:,2) N(:,:,3) N(:,:,4)];
        Me(:,:,e) = Me(:,:,e) + S_4(e,k) * R(:,:,e)' * NN(:,:,e,k)' * rho_p * NN(:,:,e,k) * R(:,:,e);
    end
     % d) Assembly to global matrices
    for j=1:6
        Idof(j,1) = 6*(Tn(e,1)-1)+j;
        Idof (6+j,1) = 6*(Tn(e,2)-1)+j;
        Idof (12+j,1) = 6*(Tn(e,3)-1)+j;
        Idof (18+j,1) = 6*(Tn(e,4)-1)+j;
    end
    K(Idof,Idof) = K(Idof,Idof) + Km(:,:,e) + Kb(:,:,e) + Ks(:,:,e);
    M(Idof,Idof) = M(Idof,Idof) + Me(:,:,e);
   
end
end