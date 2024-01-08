%=======================================
%
% Wing Modelling
%
% Master in Space and Aeronautical
% Engineering.
% Computational Enginerring CSM: Project.
% By: Jorge Simón & Iñaki Fernandez.
% Last modification: 30/12/2023
%
%=======================================
clc;
clear;
close all;
header;


%% DATA
E = 71*1e9;                           % Young modulus [Pa]
nu = 0.33;                            % Poisson's ratio
G  = E / (2 * (1 + nu));              % Shear modulus
rho = 2300;                           % Density [kg/m^3]
y_p = [0,1,0];                        % Orientation of y'
A   = 0.0062;                         % Area in [m^2]
J   = 8.949 * 1e-5;                   % Polar intertia [m^4]
Iy  = 8.801 * 1e-6;                   % Intertia Iy [m^4]
Iz  = 8.069 * 1e-5;                   % Intertia Iz [m^4]
ky  = 0.5457;                         % Shear correction factor
kz  = 0.1767;                         % Shear correction factor
kt  = 0.2920;                         % Torsion correction factor
Iyz = 0;                              % I is symmetric for y,z
shear_centre = [0.2954,0];            % In [m]
neutral_centre = [0.3083,0];          % In [m]
centre_of_mass = shear_centre;        % In [m] CHECK!
rx = 0;                               % In [m] CHECK  
ry = sqrt(Iy / A);                    % In [m] CHECK
rz = sqrt(Iz / A);                    % In [m] CHECK 
ryz = sqrt(Iyz / A);                  % In [m] CHECK
y0 = 0.3083;                          % In [m]
yc = shear_centre(1);                 % In [m]
y1 = 0.176;                           % In [m]
y2 = 0.464;                           % In [m]
h1 = 17.5*1e-3;                       % In [m] / thickness of front spar
h2 = 15*1e-3;                         % In [m] / thickness of rear spar
h3 = 6*1e-3;                          % In [m] / thickness of upper and lower surfaces

%% PREPROCESS 
disp("Loading Data For Wing")
%=======================================
%% 1.2 Mesh (discretization) data
%=======================================

% Load mesh data
load('wing.mat','xn','Tn_st','Tm_st','Tn_wb','Tm_wb','Tn_rb','Tm_rb','Tn_sk','Tm_sk','indRoot','indPointA','indPointB','indSpar1','indSpar2','n_u','n_l');
% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn_st : Nodal connectivities for beam elements (st: stringers) [Nelem x 2]
% Tn_wb, Tn_rb, Tn_sk : Nodal connectivities for shell elements (wb: wingbox, rb: ribs, sk: skin) [Nelem x 4]
% Tm_st, Tm_wb, Tm_rb, Tm_sk : Material connectivities for the different elements [Nelem x 1]
% indRoot   : Array of indices for root section nodes.
% indPointA : Index for node at point A.
% indPointB : Index for node at point B.
% indSpar1  : Array of indices for front spar centerline nodes.
% indSpar2  : Array of indices for rear spar centerline nodes.
% n_u, n_l  : Matrices containing information about the unit normals in the upper 
%             and lower surfaces, respectively. 
%       |nx1 ny1 nz1 id1|   First three columns are components of normal vector 
%       |nx2 ny2 nz2 id2|   Last (fourth) column is the nodal index
%       |      ...      |

% Number of nodes
n         = size(xn,1);                         

% Number of elements
n_elem_wb = size(Tn_wb,1);                         % Winbox
n_elem_st = size(Tn_st,1);                         % Stringers
n_elem_rb = size(Tn_rb,1);                         % Ribs
n_elem_sk = size(Tn_sk,1);                         % Skin

%=======================================
%% 1.3 Boundary conditions
%=======================================
indRoot_size = size(indRoot,1);

Up = zeros(6 * indRoot_size,3);                    % In [m]
% Define boundary conditions: Up
for i=1:6
    k = 1;
    for j=i:6:(6*indRoot_size)
        Up(j,1) = 0;
        Up(j,2) = indRoot(k);
        Up(j,3) = i;
        k = k+1;
    end
end

%=======================================
%% External forces
%=======================================

% We apply an external force in the last node in direction z.
F   = 10;                                                         % In [N]
T   = 1;                                                         % In [N]
Fa  = F * (y2-yc)/(y2-y1);                                       % In [N]
Fb  = F * (yc-y1)/(y2-y1);                                       % In [N]
FaT = -T / (y2 - y1);                                            % In [N]
FbT = T / (y2 - y1);                                             % In [N]
Fe  = [Fa indPointA 3 ; Fb indPointB 3];                       % In [N]

% We consider gravity force.
Be = zeros(n,3);
g = 0.0;                                                         % Gravity [m/s^2]
for i=1:n
    Be(i,1) = g;% Define gravity here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Be(i,2) = i;
    Be(i,3) = 3;
end

% Build pressure matrix
for i=1: size(n_u,1)
    x = xn(n_u(i,4),1);
    y = xn(n_u(i,4),2);
    pu(i) = p_u(x,y); 
end
for i=1: size(n_l,1)
    x = xn(n_l(i,4),1);
    y = xn(n_l(i,4),2);
    pl(i) = p_l(x,y); 
end

% Pressure forces
Pe = zeros(2 * size(n_u,1),3);
index = 1;
for i=1:3:(size(n_u,1) * 3)
    px = pu(index) * n_u(index,1);
    py = pu(index) * n_u(index,2);
    pz = pu(index) * n_u(index,3);
    Pe(i,:) = [-px,n_u(index,4),1];
    Pe(i+1,:) = [-py,n_u(index,4),2];
    Pe(i+2,:) = [-pz,n_u(index,4),3];
    index = index +1;
end
index = 1;
for i=(size(n_u,1)*3 + 1):3:(2 * 3 * size(n_u,1))
    px = pl(index) * n_l(index,1);
    py = pl(index) * n_l(index,2);
    pz = pl(index) * n_l(index,3);
    Pe(i,:) = [-px,n_l(index,4),1];
    Pe(i+1,:) = [-py,n_l(index,4),2];
    Pe(i+2,:) = [-pz,n_l(index,4),3];
    index = index +1;
end

% If we dont want pressure use this
Pe0 = [0,1,1];

%=======================================
%% SOLVER
%=======================================

disp("Data Loaded Succesfully!")
disp("Starting the Solver for the Wing")

%=======================================
% Part = 1 --> wingbox (shell)
% Part = 2 --> stringers (beam)
% Part = 3 --> ribs (shell)
% Part = 4 --> skin (shell)
%=======================================

% Chose part
Ndof = 6 * n;

for p=1:4
Part = p;

if Part == 1
% Compute assembly
[K_wb,M_wb,R_wb,Me_wb,NN_wb,Bs_p_wb,Bb_p_wb,Kb_wb,Ks_wb,S_4_wb,Bmn_p_wb,Bmt_p_wb] = assambly_shells(Ndof,n_elem_wb,xn,Tn_wb,Tm_wb,E,rho,h1,h2,h3,nu);

% Compute artificial rotation matrix
[K__wb,Kr_wb] = arti_rot_stiff_matr(n,n_elem_wb,Tm_wb,h1,h2,h3,xn,Tn_wb,Ndof,K_wb,E,Tn_st);

% Compute globar force matrix
[f_wb] = global_force_shell(Ndof, Fe, Be, Pe0, n_elem_wb,Tm_wb,Tn_wb,Me_wb,n,S_4_wb,R_wb,NN_wb);

%======================================================================================================================================================
elseif Part == 2

% Change values of ky,kz,kt and A
ky = 5 / 6;
kz = ky;
kt = 1;
D  = 10 * 1e-3; %[m]
A  = pi * (D / 2)^2;

% Compute assembly
[K_st,M_st,l_st,R_st,Me_st,w_st,Nm_st,Ba_p_st,Bs_p_st,Bt_p_st,Bb_p_st,Ka_st,Kb_st,Ks_st,Kt_st] = assambly_g_m(Ndof,n_elem_st,xn,Tn_st,y_p,E,A,G,rho,J,Iz,Iy,ky,kz,kt);

% Compute global force matrix
[f_st] = global_force_m(Ndof,Fe,n,Be,n_elem_st,Pe0,Me_st,w_st,l_st,R_st,Nm_st,Tn_st);

%======================================================================================================================================================
elseif Part == 3
% Change values of h
h_rb = 5 *1e-3;

% Compute assembly
[K_rb,M_rb,R_rb,Me_rb,NN_rb,Bs_p_rb,Bb_p_rb,Kb_rb,Ks_rb,S_4_rb,Bmn_p_rb,Bmt_p_rb] = assambly_shells(Ndof,n_elem_rb,xn,Tn_rb,Tm_rb,E,rho,h_rb,h_rb,h_rb,nu);

% Compute artificial rotation matrix
[K__rb,Kr_rb] = arti_rot_stiff_matr(n,n_elem_rb,Tm_rb,h_rb,h_rb,h_rb,xn,Tn_rb,Ndof,K_rb,E,Tn_st);

% Compute globar force matrix
[f_rb] = global_force_shell(Ndof, Fe, Be, Pe0, n_elem_rb,Tm_rb,Tn_rb,Me_rb,n,S_4_rb,R_rb,NN_rb);

%======================================================================================================================================================
elseif Part == 4
 
% Change values of h
h_sk = 2 * 1e-3;

% Compute assembly
[K_sk,M_sk,R_sk,Me_sk,NN_sk,Bs_p_sk,Bb_p_sk,Kb_sk,Ks_sk,S_4_sk,Bmn_p_sk,Bmt_p_sk] = assambly_shells(Ndof,n_elem_sk,xn,Tn_sk,Tm_sk,E,rho,h_sk,h_sk,h_sk,nu);

% Compute artificial rotation matrix
[K__sk,Kr_sk] = arti_rot_stiff_matr(n,n_elem_sk,Tm_sk,h_sk,h_sk,h_sk,xn,Tn_sk,Ndof,K_sk,E,Tn_st);

% Compute globar force matrix
[f_sk] = global_force_shell(Ndof, Fe, Be, Pe0, n_elem_sk,Tm_sk,Tn_sk,Me_sk,n,S_4_sk,R_sk,NN_sk);

end
end

%=======================================
%% Compute total stifness matrix
%=======================================

K_Tot = K__wb + K_st + K__rb + K__sk*1e-7;
M_Tot = M_wb + M_st + M_rb + M_sk*1e-7;

%=======================================
%% Compute total force matrix
%=======================================

f = (f_wb + f_st + f_sk*1e-7 + f_rb) ;         % Normalize

%=======================================
%% Boundary Conditions
%=======================================
% 5.1 Initialization
u = zeros(Ndof,1);

% 5.2 Prescribed and free DOFs
for p=1:size(Up)
    Ip(p) = 6*(Up(p,2)-1)+Up(p,3);        % Vcetor with prescribed degrees of freedom
    u(Ip(p),1) = Up(p,1);
end

If = setdiff(1:Ndof,Ip);

%=======================================
%% Solve system of equations(static case)
%=======================================
% 6.1 Solve system
u(If,1) = K_Tot(If,If)\(f(If,1) - K_Tot(If,Ip) * u(Ip,1));         % displacement/rotations at free DOFs
umax = max(abs(u));
fr = K_Tot * u - f;                                                % reaction forces/moments at prescribed DOFs

%=======================================
%% Computing strain and internal forces in wing elements
%=======================================

% 7.1 Get stress and strain at each Gauss point:
[sigVM_wb] = stress_VM(n_elem_wb,Tm_wb,h1,h2,h3,Tn_wb,Bb_p_wb,R_wb,u,Bmn_p_wb,Bmt_p_wb,Bs_p_wb,nu,E);
[sigVM_rb] = stress_VM(n_elem_rb,Tm_rb,h_rb,h_rb,h_rb,Tn_rb,Bb_p_rb,R_rb,u,Bmn_p_rb,Bmt_p_rb,Bs_p_rb,nu,E);
[sigVM_sk] = stress_VM(n_elem_sk,Tm_sk,h_sk,h_sk,h_sk,Tn_sk,Bb_p_sk,R_sk,u,Bmn_p_sk,Bmt_p_sk,Bs_p_sk,nu,E);

%=======================================
%% Compute Modal analysis
%=======================================

% Natural modes quantity
Nm = 6; 

% Modal analysis
[omega,phi,lambda] = mod_analy(Nm,K_Tot,M_Tot,If,Ndof);

% Obtain the frequencies
freq = omega / (2 * pi);


%========================================
%% Modal-order redcution
%========================================

% Select set of modes to project the system
Im = [1,2,3,4,5,6];

% Compute modal-order reduction
[U_s] = modal_order_reduction(Im,Ndof,f,lambda,phi);


%=======================================
%% POSTPROCESS: 2D PLOTS u
%=======================================

% Extract the values of xn(:,1) and the selected elements of u

u_z = u(3:6:end);
u_y = u(2:6:end);

% uz(1)
uz_spar1 = u_z(indSpar1);

% uz(2)
uz_spar2 = u_z(indSpar2);

% uy(1)
uy_spar1 = u_y(indSpar1);

% uy(2)
uy_spar2 = u_y(indSpar2);

% xn values
xn_values = xn(indSpar1);

% thetax average
thetax = (uz_spar2 - uz_spar1)/(y2-y1);

% uz average
uz_average = uz_spar1 + thetax*(yc-y1);

% uy average
uy_average = (uy_spar1 + uy_spar2) / 2;

% Plotting parameters
fs = 20;

% Plotting the values
figure;
plot(xn_values, round(uz_average,5,"decimals")*1e3, '-', 'LineWidth', 2, Color='red');
%plot(xn_values, uz_average*1e3, '-', 'LineWidth', 2, Color='red');

% Adding labels and title
xlabel('Position [m]', 'Interpreter','latex',FontSize=fs);
ylabel('Displacement $u_{z}$ [mm]','Interpreter','latex',FontSize=fs);
title('Position VS Displacement', 'Interpreter','latex',FontSize=fs);
set(gca,"ticklabelinterpreter",'latex')

% Display the grid
grid on;

% Plotting the values
figure;
%plot(xn_values, round(thetax,5,"decimals"), '-', 'LineWidth', 2, Color='red');
plot(xn_values, thetax, '-', 'LineWidth', 2, Color='red');

% Adding labels and title
xlabel('Position [m]', 'Interpreter','latex',FontSize=fs);
ylabel('Displacement $\theta_{x}$ [Rad]','Interpreter','latex',FontSize=fs);
title('Position VS Displacement', 'Interpreter','latex',FontSize=fs);
set(gca,"ticklabelinterpreter",'latex')

% Display the grid
grid on;

%=======================================
%% POSTPROCESS: 2D PLOTS u*
%=======================================

% Extract the values of xn(:,1) and the selected elements of u

u_z = U_s(3:6:end);
u_y = U_s(2:6:end);

% uz(1)
uz_spar1 = u_z(indSpar1);

% uz(2)
uz_spar2 = u_z(indSpar2);

% uy(1)
uy_spar1 = u_y(indSpar1);

% uy(2)
uy_spar2 = u_y(indSpar2);

% xn values
xn_values = xn(indSpar1);

% thetax average
thetax = (uz_spar2 - uz_spar1)/(y2-y1);

% uz average
uz_average = uz_spar1 + thetax*(yc-y1);

% uy average
uy_average = (uy_spar1 + uy_spar2) / 2;

% Plotting parameters
fs = 20;

% Plotting the values
figure;
plot(xn_values, uz_average*1e3, '-', 'LineWidth', 2, Color='red');

% Adding labels and title
xlabel('Position [m]', 'Interpreter','latex',FontSize=fs);
ylabel('Displacement $u^{*}_{z}$ [mm]','Interpreter','latex',FontSize=fs);
title('Position VS Displacement', 'Interpreter','latex',FontSize=fs);
set(gca,"ticklabelinterpreter",'latex')

% Display the grid
grid on;

% Plotting the values
figure;
plot(xn_values, thetax, '-', 'LineWidth', 2, Color='red');

% Adding labels and title
xlabel('Position [m]', 'Interpreter','latex',FontSize=fs);
ylabel('Displacement $\theta^{*}_{x}$ [Rad]','Interpreter','latex',FontSize=fs);
title('Position VS Displacement', 'Interpreter','latex',FontSize=fs);
set(gca,"ticklabelinterpreter",'latex')

% Display the grid
grid on;

%=======================================
%% POSTPROCESS: 3D PLOTS
%======================================= 

% Additional plot functions useful to visualize 3D model and modes

scale = 10; 
plotDeformed('wing',xn,Tn_wb,u,scale,sigVM_wb); % For wingbox elements
plotDeformed('wing',xn,Tn_rb,u,scale,sigVM_rb); % For rib elements
plotDeformed('wing',xn,Tn_sk,u,scale,sigVM_sk); % For skin elements
% % This function plots the deformed structure and Von Mises stress distribution: 
% % xn : Nodal coordinates matrix [Nnodes x 3]
% % Tn : Nodal connectivities matrix [Nelem x 4]
% % u : Displacements vector obtained from the system solution. It is expected
% %     to be given as a column vector with dimensions [Ndof x 1].
% % scale : Scale factor to amplify the displacements (set to appropriate 
% %         number to visualize the deformed structure properly).
% % sigVM : Von Mises stress at each Gauss point. It is expected to be given as 
% %         a matrix with dimensions [Nelem x Ngauss].
% 
imodes = [1,2,3,4,5,6];
plotModes('wing',phi,freq,imodes)
% % This function plots the specified modes resulting from a modal analysis
% % in sets of 9.
% % Phi : Modal displacements matrix in which each column corresponds to the
% %       mode shape of the corresponding mode. Expected dimensions [Ndof x Nmodes]
% % freq : Natural frequencies array. Expected dimensions [Nmodes x 1]
% % imodes : Array selecting which modes to plot. A maximum of 9 modes must
% %          can be selected. Example: imodes = [1,2,3,4,5,6,7,8,9] will plot
% %          the modes stored in the first 9 columns of Phi / imodes = [1,4,5,10] 
% %          will plot modes in columns 1, 4, 5 and 10 of Phi. 