%=======================================
%
% Shell Modelling
%
% Master in Space and Aeronautical
% Engineering.
% Computational Enginerring CSM: Project.
% By: Jorge Simón & Iñaki Fernandez.
% Last modification: 20/12/2023
%
%=======================================
clc;
clear;
close all;
format long;
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
disp("Loading Data For Shells")
%=======================================
%% 1.2 Mesh (discretization) data
%=======================================

% Load mesh data
load('shell.mat','xn','Tn','Tm','indRoot','indPointA','indPointB','indSpar1','indSpar2');
% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn : Nodal connectivities [Nelem x 4]
% Tm : Material connectivities [Nelem x 1]
% indRoot   : Array of indices for root section nodes.
% indPointA : Index for node at point A.
% indPointB : Index for node at point B.
% indSpar1  : Array of indices for front spar centerline nodes.
% indSpar2  : Array of indices for rear spar centerline nodes.

n_elem = size(Tn,1);                         % Number of elements
n      = size(xn,1);                         % Number of nodes 
% Define boundary conditions and forces data matrices: Up, Fe, Pe, Be
% ...
%=======================================
%% 1.3 Boundary conditions
%=======================================

Up = zeros(6*28,3);                    % In [m]
% Define boundary conditions: Up
for i=1:6
    k = 1;
    for j=i:6:(6*28)
        Up(j,1) = 0;
        Up(j,2) = indRoot(k);
        Up(j,3) = i;
        k = k+1;
    end
end

%=======================================
%% 1.4 External forces
%=======================================
% We apply an external force in the last node in direction z.
F   = 1;                                                         % In [N]
T   = 1;                                                         % In [N]
Fa  = F * (y2-yc)/(y2-y1);                                       % In [N]
Fb  = F * (yc-y1)/(y2-y1);                                       % In [N]
FaT = -T / (y2 - y1);                                            % In [N]
FbT = T / (y2 - y1);                                             % In [N]
Fe  = [FaT indPointA 3 ; FbT indPointB 3];                       % In [N]

% We consider gravity force.
Be = zeros(n,3);
g = 0;                    % Gravity [m/s^2]
for i=1:n
    Be(i,1) = g;               % Define gravity here
    Be(i,2) = i;
    Be(i,3) = 3;
end

% Pressure forces
Pe=zeros(3,3);

%% SOLVER
disp("Data Loaded Succesfully!")
disp("Starting the Solver for the Shells")

%=======================================
%% Assembly of global matrices
%=======================================
% Define degrees of freedom
Ndof = 6*n;

% Obtain system matrices
[K,M,R,Me,NN,Bs_p,Bb_p,Kb,Ks,S_4,Bmn_p,Bmt_p] = assambly_shells(Ndof,n_elem,xn,Tn,Tm,E,rho,h1,h2,h3,nu);

%=======================================
%% Compute artificial rotation stiffness matrix
%=======================================

[K_,Kr,s,S,k__p,n_mat] = artificial_stiffness(K,Ndof,n_elem,n,Tn,xn,Tm,E,h1,h2,h3);

%=======================================
%% Compute global force vector
%=======================================

[P,B,f,be,pe,fe] = global_force(Fe,Be,n_elem,n,Ndof,h1,h2,h3,Tn,Tm,Me,S_4,R,NN);

%=======================================
%% Boundary Conditions
%=======================================

% 5.1 Initialization
u = zeros(Ndof,1);

% 5.2 Prescribed and free DOFs
for p=1:size(Up)
    Ip(p) = 6*(Up(p,2)-1)+Up(p,3);                   % Vector with prescribed degrees of freedom
    u(Ip(p),1) = Up(p,1);
end

If = setdiff(1:Ndof,Ip);

%=======================================
%% Solve system equations (static case)
%=======================================

% 6.1 Solve system
u(If,1) = K_(If,If)\(f(If,1) - K_(If,Ip) * u(Ip,1)); % displacement/rotations at free DOFs
umax = max(abs(u));
fr = K_ * u - f;                                     % reaction forces/moments at prescribed DOFs

%=======================================
%% Computing local strain and stress in shell elements
%=======================================

[Idof,Eb_p,Em_p,Es_p,Cp,Cs,sig_m,sig_s,sig_b,sig_VM] = postprocess_shells(n_elem,h1,h2,h3,Tm,Tn,Bb_p,Bmn_p,Bmt_p,Bs_p,R,u,nu,E);

%=======================================
%% Compute Modal analysis
%=======================================

Nm = 6; % Natural modes quantity

% Modal analysis
[omega,phi,lambda] = mod_analy(Nm,K_,M,If,Ndof);

% Obtain frequencies
frec = omega / (2*pi);

% Obtain vibration modes for uz,uy and thetax for each DOF
uz_modes     = zeros(size(xn,1),6);
uy_modes     = zeros(size(xn,1),6);
thetax_modes = zeros(size(xn,1),6);

for i=1:6
    uz_modes(:,i)     = phi(3:6:end,i);
    uy_modes(:,i)     = phi(2:6:end,i);
    thetax_modes(:,i) = phi(4:6:end,i);
end
for i=1:length(indSpar1)
    for j=1:6
        thetax_av_modes(i,j) = (uz_modes(indSpar2(i),j) - uz_modes(indSpar1(i),j))/(y2 - y1);
        uz_av_modes(i,j) = uz_modes(indSpar1(i),j) + thetax_modes(i,j) * (yc - y1);
        uy_av_modes(i,j) = (uy_modes(indSpar2(i),j) + uy_modes(indSpar1(i),j)/2);
    end
end

disp("Solver Finished Succesfully!")
%=======================================
%% POSTPROCESS: SAVE DATA AND PLOTS 
%=======================================
disp("Starting Postporcess")

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
ylabel('Displacement $\theta{x}$ [Rad]','Interpreter','latex',FontSize=fs);
title('Position VS Displacement', 'Interpreter','latex',FontSize=fs);
set(gca,"ticklabelinterpreter",'latex')

% Display the grid
grid on;

% % Specify the folder name
% folderName = 'modes_plots';
% 
% % Create the folder if it doesn't exist
% if ~exist(folderName, 'dir')
%     mkdir(folderName);
% end
% 
% for i = 1:size(uz_modes,2)
%     % Create a new figure for each iteration
%     figure;
% 
%     % Plotting the values
%     plot(xn_values, round(uz_av_modes(:,i) * 1e3,4,"decimals"), '-', 'LineWidth', 2, 'Color', 'blue');
% 
%     % Adding labels and title
%     xlabel('Position [m]', 'Interpreter', 'latex', 'FontSize', fs);
%     ylabel('Displacement $u_{z}$ [mm]', 'Interpreter', 'latex', 'FontSize', fs);
%     title(['Mode ' num2str(i)], 'Interpreter', 'latex', 'FontSize', fs);
%     set(gca, 'ticklabelinterpreter', 'latex');
% 
%     % Display the grid
%     grid on;
% 
%     % Save the figure in the specified folder
%     fileName = fullfile(folderName, ['uz_Mode_' num2str(i) '_Plot.png']);
%     saveas(gcf, fileName);
% 
% end
% 
% for i = 1:size(uy_modes,2)
%     % Create a new figure for each iteration
%     figure;
% 
%     % Plotting the values
%     plot(xn_values, round(uy_av_modes(:,i) * 1e3,4,"decimals"), '-', 'LineWidth', 2, 'Color', 'green');
% 
%     % Adding labels and title
%     xlabel('Position [m]', 'Interpreter', 'latex', 'FontSize', fs);
%     ylabel('Displacement $u_{y}$ [mm]', 'Interpreter', 'latex', 'FontSize', fs);
%     title(['Mode ' num2str(i)], 'Interpreter', 'latex', 'FontSize', fs);
%     set(gca, 'ticklabelinterpreter', 'latex');
% 
%     % Display the grid
%     grid on;
% 
%     % Save the figure in the specified folder
%     fileName = fullfile(folderName, ['uy_Mode_' num2str(i) '_Plot.png']);
%     saveas(gcf, fileName);
% 
% end
% 
% for i = 1:size(thetax_modes,2)
%     % Create a new figure for each iteration
%     figure;
% 
%     % Plotting the values
%     plot(xn_values, round(thetax_av_modes(:,i),4,"decimals"), '-', 'LineWidth', 2, 'Color', 'magenta');
% 
%     % Adding labels and title
%     xlabel('Position [m]', 'Interpreter', 'latex', 'FontSize', fs);
%     ylabel('Displacement $\theta_{x}$ [Rad]', 'Interpreter', 'latex', 'FontSize', fs);
%     title(['Mode ' num2str(i)], 'Interpreter', 'latex', 'FontSize', fs);
%     set(gca, 'ticklabelinterpreter', 'latex');
% 
%     % Display the grid
%     grid on;
% 
%     % Save the figure in the specified folder
%     fileName = fullfile(folderName, ['thetax_Mode_' num2str(i) '_Plot.png']);
%     saveas(gcf, fileName);
% 
% end

disp("Postprocess Finished Succesfully!")