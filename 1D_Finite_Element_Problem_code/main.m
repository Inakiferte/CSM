%=======================================
% Master in Space and Aeronautical
% Engineering.
% Computational Engineering: Assignment 1
% CSM.
% By: Jorge Simón & Iñaki Fernandez Tena
% Last modification: 26/11/2023.
%=======================================
clear; 
close all;
clc

% !!!!!!!!!!!!!!! PLEASE RUN ALL THE FUNCTIONS BEFORE RUNING THE main.m CODE !!!!!!!!!!!!!!!!!!
addpath = "C:\Users\User\Desktop\Iñaki\Fisika\Master\MASE\Courses\Computational_Engineering\Assigments\1D_FEM_Elements\code\functions\";
%=======================================
% Define input data
%=======================================

% Physical input data
A = 120e-6;                              % m^2 Area
E = 75e9;                                % Pa young modulus
L = 2;                                   % m Lenght
c = 1e4;                                 % Nm^-1
lambda = 2;                              % Constant

% Numerical input data
h    = [L / 2, L / 4, L / 8, L / 16, L / 32, L / 64, L / 128];

for a=2                                % Change a = 2 or 3 for linear and quadratic solver respectively
for p=1: length(h)
n_el = L / h(p);                       % Total number of elements
n_ne = a;                              % Number of nodes per element
n    = n_ne + (n_el - 1)*(n_ne - 1) ;  % Total number of nodes
n_G  = a;                              % Model (2 = linear; 3 = quadratic)


% Define matrixes

x  = zeros(n,1);                        % Nodal coordinates array
Tn = zeros(n_el,n_ne);                  % Nodal cennectivity table

% Fill x and Tn
if a == 2
    for i=2: n
        x(1) = 0;
        x(i) = x(i - 1) + h(p);
    end
else 
    for i=2: n
    x(1) = 0;
    x(i) = x(i - 1) + h(p)/2;
    end
end
for e=1: n_el
    for i=1: n_ne
        Tn(e,i) = (e - 1)*(n_ne - 1) + i;
    end
end

%=======================================
%% 2- Computation of the global stiffness
%=======================================

B  = zeros(1,n_ne,n_el,n_G);
Ke = zeros(n_ne,n_ne,n_el);
K  = zeros(n,n);
dN_dxhi = zeros(1,n_G);

for e=1: n_el
    % a) Get element's size
    l_e = x(Tn(e, n_ne)) - x(Tn(e, 1));
    xp1  = x(Tn(e, 1));
    xp2  = x(Tn(e, 2));
    xp3  = x(Tn(e, n_ne));
    for k=1:n_G % (see appendix B)
        [wk,xhi] = gauss(k, n_G);
        % b) Compute B matrix evaluated at xhi_k (see Appendix A)
        [J, N, dN_dxhi] = shape(e,xhi,n_G,xp1,xp2,xp3);
        for i=1: n_ne
            B(1,i,e,k) = (1 / J) * dN_dxhi(i);
        end
        % c) Compute the element stiffness matrix evaluated at xhi_k
        Ke(:,:,e) = Ke(:,:,e) + wk * J * E * A * (transpose(B(:,:,e,k)) * B(:,:,e,k));
    end
    % d) Assembly to global stiffness matrix
    for i=1: n_ne
        for j=1: n_ne
            K(Tn(e,i),Tn(e,j)) = K(Tn(e,i),Tn(e,j)) + Ke(i,j,e);
        end
    end
end

%=======================================
%% 3- Computation of the global force vector
%=======================================

N_matrix   = zeros(1,n_ne,n_el,n_G);
x_G = zeros(n_G,n_el);
fe  = zeros(n_ne,n_el);
f   = zeros(n,1);

for e=1: n_el
    % a) Get nodal coordinates and element%s size
    xe  = x(Tn(e,:));
    l_e = x(Tn(e, n_ne)) - x(Tn(e, 1));
    xp1  = x(Tn(e, 1));
    xp2  = x(Tn(e, 2));
    xp3  = x(Tn(e, n_ne));
    for k=1: n_G
        [wk,xhi] = gauss(k, n_G);
        % b) Compute 𝐍𝐍 matrix evaluated at 𝜉𝑘 (see Appendix A)
        [J, N, dN_dxhi] = shape(e,xhi,n_G,xp1,xp2,xp3);
        for i=1: n_ne
            N_matrix(1,i,e,k) = N(i);
        end
        % c) Compute Gauss point coordinate and element force vector evaluated at 𝜉𝑘
        x_G(k,e) = (N_matrix(:,:,e,k)) * xe;
        % d) Evaluate body force at Gauss point
        bk = body_force(x_G(k,e),c,lambda);
        %e) Compute element force vector evaluated at 𝜉𝑘
        fe(:,e) = fe(:,e) + wk * J * ((N_matrix(:,:,e,k))') * bk;
    end
    % d) Assembly to global force matrix
    for i=1: n_ne
        f(Tn(e,i),1) = f(Tn(e,i),1) + fe(i,e);
    end
end
%=======================================
%% 4- Global system of equations
%=======================================

% a) Apply conditions.
nu_r = [1];
nu_f = setdiff(1:n,nu_r);

% b) Partitioned system of equations

K_ff = K(nu_f,nu_f);
K_fr = K(nu_f,nu_r);
K_rf = K(nu_r,nu_f);
K_rr = K(nu_r,nu_r);

% c) System resolution.
ur        = [0];
u(nu_r,1) = ur;
u(nu_f,1) = (inv(K_ff)) * (f(nu_f,1) - K_fr * u(nu_r,1));
fr        = K_rr * u(nu_r,1) + K_rf * u(nu_f,1) - f(nu_r,1);

%=======================================
%% 5- Computation of stress
%=======================================

sigG = zeros(n_G,n_el);
if a==2
    for e=1: n_el
        % a) Obtain element nodes displacement
        ue_linear(:,e) = u(Tn(e,:));
        % b) Obtain stress at Gauss point
        for k=1: n_G
            sigG(k,e) = E * (B(:,:,e,k)) * ue_linear(:,e);
        end
    end
else
    for e=1: n_el
        % a) Obtain element nodes displacement
        ue_quadratic(:,e) = u(Tn(e,:));
        % b) Obtain stress at Gauss point
        for k=1: n_G
            sigG(k,e) = E * (B(:,:,e,k)) * ue_quadratic(:,e);
        end
    end

end

%=======================================
%% 7- Save info
%=======================================
if p==1
    u1 = u;
    x1 = x;
    sig1 = interp1(x_G(:), sigG(:),x,"linear","extrap");
end

if p==2
    u2 = u;
    x2 = x;
    sig2 = interp1(x_G(:), sigG(:),x,"linear","extrap");
end

if p==3
    u3 = u;
    x3 = x;
    sig3 = interp1(x_G(:), sigG(:),x,"linear","extrap");
end

if p==4
    u4 = u;
    x4 = x;
    sig4 = interp1(x_G(:), sigG(:),x,"linear","extrap");
end

if p==5
    u5 = u;
    x5 = x;
    sig5 = interp1(x_G(:), sigG(:),x,"linear","extrap");
end

if p==6
    u6 = u;
    x6 = x;
    sig6 = interp1(x_G(:), sigG(:),x,"linear","extrap");
end

if p==7
    u7 = u;
    x7 = x;
    sig7 = interp1(x_G(:), sigG(:),x,"linear","extrap");
end
 end

%=======================================
%% Compute the analytical function
%=======================================
x_analy = linspace(0,2,1000);
u_analy = u_analytical(x_analy,c,E,A,lambda,L);
sigma_analy = sigma_analytical(x_analy,c,A,lambda,L);

%=======================================
%% 8- Compute the error
%=======================================
e_u_linear = zeros(length(h),1);
e_s_linear = zeros(length(h),1);
e_u_quadratic = zeros(length(h),1);
e_s_quadratic = zeros(length(h),1);
if a==2
for p=1:length(h)
    if p == 1
    [e_u_linear(p), e_s_linear(p)] = error(u1, u_analy, sig1, sigma_analy);
    end
    if p == 2
    [e_u_linear(p), e_s_linear(p)] = error(u2, u_analy, sig2, sigma_analy);
    end
    if p == 3
    [e_u_linear(p), e_s_linear(p)] = error(u3, u_analy, sig3, sigma_analy);
    end
    if p == 4
    [e_u_linear(p), e_s_linear(p)] = error(u4, u_analy, sig4, sigma_analy);
    end
    if p == 5
    [e_u_linear(p), e_s_linear(p)] = error(u5, u_analy, sig5, sigma_analy);
    end
    if p == 6
    [e_u_linear(p), e_s_linear(p)] = error(u6, u_analy, sig6, sigma_analy);
    end
    if p == 7
    [e_u_linear(p), e_s_linear(p)] = error(u7, u_analy, sig7, sigma_analy);
    end
end
end
if a ==3
for p=1:length(h)
    if p == 1
    [e_u_quadratic(p), e_s_quadratic(p)] = error(u1, u_analy, sig1, sigma_analy);
    end
    if p == 2
    [e_u_quadratic(p), e_s_quadratic(p)] = error(u2, u_analy, sig2, sigma_analy);
    end
    if p == 3
    [e_u_quadratic(p), e_s_quadratic(p)] = error(u3, u_analy, sig3, sigma_analy);
    end
    if p == 4
    [e_u_quadratic(p), e_s_quadratic(p)] = error(u4, u_analy, sig4, sigma_analy);
    end
    if p == 5
    [e_u_quadratic(p), e_s_quadratic(p)] = error(u5, u_analy, sig5, sigma_analy);
    end
    if p == 6
    [e_u_quadratic(p), e_s_quadratic(p)] = error(u6, u_analy, sig6, sigma_analy);
    end
    if p == 7
    [e_u_quadratic(p), e_s_quadratic(p)] = error(u7, u_analy, sig7, sigma_analy);
    end
end
end
end
%=======================================
%% 9- Plot stress
%=======================================


%% Create the plots
figure(1);
fs = 25;

plot(x1, u1, LineWidth=2, Color='blue');
hold on
plot(x2, u2, LineWidth=2, Color='green');
plot(x3, u3, LineWidth=2, Color='red');
plot(x4, u4, LineWidth=2, Color='magenta');
plot(x5, u5, LineWidth=2, Color='yellow');
plot(x6, u6, LineWidth=2, Color='cyan');
plot(x7, u7, LineWidth=2, Color='black');
plot(x_analy,u_analy,'--',LineWidth=2, Color='red')
xlabel('Bar position(m)','Interpreter','latex',FontSize=fs)
ylabel(' Displacement (m)','Interpreter','latex',FontSize=fs)
title('Bar Displacement vs Distance','Interpreter','latex', FontSize=fs)
legend('h=L/2', 'h=L/4', 'h=L/8', 'h=L/16', 'h=L/32', 'h=L/64', 'h=L/128', 'Analytical', fontsize=fs-15)
set(legend, 'Interpreter','latex', 'Location', 'North West'); 
grid;
%axis([0 2 0 0.05])

figure(2);
plot(x1, sig1*1e-9, LineWidth=2, Color='blue');
hold on
plot(x2, sig2*1e-9, LineWidth=2, Color='green');
plot(x3, sig3*1e-9, LineWidth=2, Color='red');
plot(x4, sig4*1e-9, LineWidth=2, Color='magenta');
plot(x5, sig5*1e-9, LineWidth=2, Color='yellow');
plot(x6, sig6*1e-9, LineWidth=2, Color='cyan');
plot(x7, sig7*1e-9, LineWidth=2, Color='black');
plot(x_analy,sigma_analy*1e-9,'--',LineWidth=2, Color='red');
xlabel('Bar position(m)','Interpreter','latex',FontSize=fs)
ylabel(' Stress (GPa)','Interpreter','latex',FontSize=fs)
title('Bar Stress vs Distance','Interpreter','latex', FontSize=fs)
legend('h=L/2', 'h=L/4', 'h=L/8', 'h=L/16', 'h=L/32', 'h=L/64', 'h=L/128', 'Analytical', fontsize=fs-15)
set(legend, 'Interpreter','latex', 'Location', 'North West'); 
grid;

figure(3)
plot(log(h), log(e_u_linear), LineWidth=2, Color='blue');
hold on
plot(log(h), log(e_u_quadratic), LineWidth=2, Color='green');
xlim([log(h(end)), log(h(1))]);
xlabel('Log(h)','Interpreter','latex',FontSize=fs)
ylabel('Log(e)','Interpreter','latex',FontSize=fs)
title('Linear','Interpreter','latex', FontSize=fs)
grid;

disp( "=======================================================");
disp( "1D Finite Element Problem.");
disp( "Master in Space and Aeronautical Engineering.");
disp( "Computational Engineering Task 01.");
disp( "By: Jorge Simón & Iñaki Fernandez.");
disp( "=======================================================");
disp( " ");
disp( "=======================================================");
disp( "The program has finished successfully!");
disp( "You can find the result plotted in your screen for: Error, Stress and Displacement.");
disp( "Results from linear to quadratic solver can be changed in line 28 from a = 2 or 3 for linear and quadratic solver respectively.")
disp( "=================================================");
