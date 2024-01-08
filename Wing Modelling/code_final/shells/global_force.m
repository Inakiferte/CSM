function [P,B,f,be,pe,fe] = global_force(Fe,Be,n_elem,n,Ndof,h1,h2,h3,Tn,Tm,Me,S_4,R,NN)
%% 4) Compute global force vector
% 4.1) Point loads
f = zeros(Ndof,1);
for q=1:size(Fe)
    f(6*(Fe(q,2)-1)+Fe(q,3),1) = f(6*(Fe(q,2)-1)+Fe(q,3),1) + Fe(q,1); 
end

% 4.2 Nodal distributed forces
P = zeros(n,6);

%for r=1:size(Pe)
%    P(Pe(r,2),Pe(r,3)) = P(Pe(r,2),Pe(r,3)) + Pe(r,1);  % Da fallo porque Pe esta como matriz de ceros, hay que rellenarla, mirar Dist_F_Det (Xavi)
%end

% 4.3 Nodal body forces
B = zeros (n,6);
for s=1:size(Be)
    B(Be(s,2),Be(s,3)) = B(Be(s,2),Be(s,3)) + Be(s,1);
end

% 4.4 Assembly process
%fe = zeros(6*4,Ne_s);
for e=1:n_elem
    if Tm(e) ==1
        h = h1;
    elseif Tm(e) == 2
        h = h2;
    else 
        h = h3;
    end
    % a) Compute element force vector
    be(:,e) = [B(Tn(e,1),:), B(Tn(e,2),:), B(Tn(e,3),:), B(Tn(e,4),:)]';
    pe(:,e) = [P(Tn(e,1),:), P(Tn(e,2),:), P(Tn(e,3),:), P(Tn(e,4),:)]';
    fe(:,e) = Me(:,:,e) * be(:,e);
    for k=1:4
        fe(:,e) = fe(:,e) + S_4(e,k) * R(:,:,e)' * NN(:,:,e,k)' * NN(:,:,e,k) * R(:,:,e) * pe(:,e);
    end

    % b) Assembly to global force vector
    for j=1:6
        Idof(j,1) = 6*(Tn(e,1)-1)+j;
        Idof(j+6,1) = 6*(Tn(e,2)-1)+j;
        Idof(j+12,1) = 6*(Tn(e,3)-1)+j;
        Idof(j+18,1) = 6*(Tn(e,4)-1)+j;
    end
    f(Idof,1) = f(Idof,1) + fe(:,e);
end
end