function [K_,Kr,s,S,k__p,n_mat] = artificial_stiffness(K,Ndof,n_elem,n,Tn,xn,Tm,E,h1,h2,h3)
%% 3) Compute artificial rotation stiffness matrix
% 3.1 Find nodal normal to set criteria for finding coplanar nodes
n_mat = zeros (3,n);
for e=1:n_elem
    if Tm(e) ==1
        h = h1;
    elseif Tm(e) == 2
        h = h2;
    else 
        h = h3;
    end
    % a) Compute normal and surface
    s = cross((xn(Tn(e,3),:)'- xn(Tn(e,1),:)'),(xn(Tn(e,4),:)'- xn(Tn(e,2),:)'))/2;
    S(e) = sqrt(s(1)^2 + s(2)^2 + s(3)^2);
    k__p(:,e) = s/S(e);

    % b) Assemble to get nodal normal
    for i=1:4 
        n_mat(:,Tn(e,i)) = n_mat(:,Tn(e,i)) + k__p(:,e);
    end
end

% 3.2 Compute artificial rotation matrix
Kr = sparse(Ndof,Ndof);
for e=1:n_elem
    if Tm(e) ==1
        h = h1;
    elseif Tm(e) == 2
        h = h2;
    else 
        h = h3;
    end
    for i=1:4
        alpha = acos( dot(n_mat(:,Tn(e,i)), k__p(:,e)) / norm(n_mat(:,Tn(e,i))));
        if rad2deg(alpha) < 5              % ¡¡¡¡REVISAR!!!!
            Idof = 6*(Tn(e,i)-1) + [4; 5; 6];
            Kr(Idof,Idof) = Kr(Idof,Idof) + E * h * S(e) * k__p(:,e) * k__p(:,e)';     
        end
    end
end

% 3.3 Update stiffness matrix
K_ = sparse(Ndof,Ndof);
K_ = K + Kr;
end