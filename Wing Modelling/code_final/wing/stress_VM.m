function [sig_VM] = stress_VM(n_elem,Tm,h1,h2,h3,Tn,Bb_p,R,u,Bmn_p,Bmt_p,Bs_p,nu,E)
%% 7) Postprocess: Computing local strain and stress in shell elements
% 7.1 Get stress and strain at each Gauss point:
for e=1:n_elem
    if Tm(e) ==1
        h = h1;
    elseif Tm(e) == 2
        h = h2;
    else 
        h = h3;
    end
    % a) Get each strain component:
    for j=1:6
        Idof(j,1) = 6*(Tn(e,1)-1)+j;
        Idof(j+6,1) = 6*(Tn(e,2)-1)+j;
        Idof(j+12,1) = 6*(Tn(e,3)-1)+j;
        Idof(j+18,1) = 6*(Tn(e,4)-1)+j;
    end
    for k=1:4
        Eb_p(:,e,k) = Bb_p(:,:,e,k) * R(:,:,e) * u(Idof,1);
        Em_p(1:2,e,k) = Bmn_p(:,:,e,k) * R(:,:,e) * u(Idof,1);
        Em_p(3,e,k) = Bmt_p(:,:,e) * R(:,:,e) * u(Idof,1);
        Es_p(:,e,k) = Bs_p(:,:,e) * R(:,:,e) * u(Idof,1);
    end
    
    % b) Get stress:
    Cp = [1 nu 0; nu 1 0; 0 0 (1 - nu)/2]*E/(1-(nu)^2);
    Cs = [1 0; 0 1]*E/(2*(1+nu));
    for k=1:4
        sig_m(:,e,k) = Cp * Em_p(:,e,k);
        sig_s(:,e,k) = Cs * Es_p(:,e,k);
        sig_b(:,e,k) = Cp * h * Eb_p(:,e,k)/2;

        sig_mas = [sig_m(:,e,k) + sig_b(:,e,k); sig_s(:,e,k)]';
        sig_vm_mas = sqrt(sig_mas(1)^2 + sig_mas(2)^2 - sig_mas(1)*sig_mas(2) + 3*(sig_mas(3)+abs(sig_mas(4))+abs(sig_mas(5))));

        sig_menos = [sig_m(:,e,k) - sig_b(:,e,k); sig_s(:,e,k)]';
        sig_vm_menos = sqrt(sig_menos(1)^2 + sig_menos(2)^2 - sig_menos(1)*sig_menos(2) + 3*(sig_menos(3)+abs(sig_menos(4))+abs(sig_menos(5))));
        
        sig_VM (e,k) = max ([sig_vm_menos sig_vm_mas]);
    end
end
end