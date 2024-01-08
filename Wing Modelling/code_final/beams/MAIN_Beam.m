%=======================================
%
% Beam Modelling
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
centre_of_mass = shear_centre;        % In [m] 
rx = 0;                               % In [m]   
ry = sqrt(Iy / A);                    % In [m] 
rz = sqrt(Iz / A);                    % In [m]  
ryz = sqrt(Iyz / A);                  % In [m] 


%% PREPROCESS 
disp("Loading Data For Beams")
%=======================================
%% 1.2 Mesh (discretization) data
%=======================================

% Load mesh data
load('beam.mat','xn','Tn','Tm');

% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn : Nodal connectivities [Nelem x 2]
% Tm : Material connectivities [Nelem x 1]
n_elem = size(Tn,1);                         % Number of elements
n      = size(xn,1);                         % Number of nodes 

%=======================================
%% 1.3 Boundary conditions
%=======================================

% Define boundary conditions: Up
Up = [0 1 1; 0 1 2; 0 1 3; 0 1 4; 0 1 5; 0 1 6]; % In [m]

%=======================================
%% 1.4 External forces
%=======================================
% We apply an external forces in the last node in direction z.
Fe = [1,n,4];                                    % In [N]

% There is no distribution of the load since force is pointlike.
Qe = [0,0,0];

% We consider no gravity force.
Be = zeros(n,3);
for i=1:n
    Be(i,1) = 0.0;             % Define gravity here
    Be(i,2) = i;
    Be(i,3) = 3;
end

%% SOLVER
disp("Data Loaded Succesfully!")
disp("Starting the Solver for the Beam")
%=======================================
%% 2- Assembly of global matrices
%=======================================

% 2.1- Initialization
N_dof = 6 * n;                          % Total degrees of freedom
K     = sparse(N_dof,N_dof);
M     = sparse(N_dof,N_dof);            % Porposed in modal analysis instead of zeros

% To avoid recomputing the system matrices use a save/load structure:
if 1 % Set to 1 to (re)compute the system matrices and '0' to load them
    % 2.2 Assembly process
    [K,M,l,R,Me,w,Nm,Ba_p,Bs_p,Bt_p,Bb_p,Ka,Kb,Ks,Kt] = assambly_g_m(N_dof,n_elem,xn,Tn,y_p,E,A,G,rho,J,Iz,Iy,ky,kz,kt);
    save('beam_matrices.mat','K','M');
    
else
    
    % Load previously computed results
    load('beam_matrices.mat','K','M');
    
end
%=======================================
%% 3- Compute global force vector
%=======================================

[f_hat] = global_force_m(N_dof,Fe,n,Be,n_elem,Qe,Me,w,l,R,Nm,Tn);

%=======================================
%% 4- Boundary Conditions
%=======================================
% 4.1 Initialization
u_hat = zeros(N_dof,1);

% 4.2 Prescribed and free DOFs
for p=1:size(Up,1)
    I_p(p) = 6 * (Up(p,2) - 1) + Up(p,3);
    u_hat(I_p(p),1) = Up(p,1);
end
If = setdiff(1:N_dof,I_p);

%=======================================
%% 5- Solve system of equations(static case)
%=======================================
% 5.1 Solve system
u_hat(If,1) = inv(K(If,If)) * (f_hat(If,1) - K(If,I_p)*u_hat(I_p,1));
f_hatR = K * u_hat - f_hat;

%=======================================
%% 6- Computing strain and internal forces in beam elements
%=======================================

[epsilon_a,epsilon_s,epsilon_t,epsilon_b,Fx,Fy,Fz,Mx,My,Mz] = strain_and_forces(n_elem,Tn,u_hat,R,Ba_p,Bs_p,Bt_p,Bb_p,Ka,Kb,Ks,Kt);

%=======================================
%% 7- Compute Modal analysis
%=======================================
Nm = 6;                                 % Natural modes quantity

% Modal analysis
[omega,phi] = mod_analy(Nm,K,M,If,N_dof);

% Obtain the frequencies
freq        = omega / (2 * pi);

% Obtain vibration modes for uz,uy and thetax for each DOF
uz_modes     = zeros(size(xn,1),6);
uy_modes     = zeros(size(xn,1),6);
thetax_modes = zeros(size(xn,1),6);

for i=1:6
    uz_modes(:,i)     = phi(3:6:end,i);
    uy_modes(:,i)     = phi(2:6:end,i);
    thetax_modes(:,i) = phi(4:6:end,i);
end

disp("Solver Finished Succesfully!")
%=======================================
%% POSTPROCESS: SAVE DATA AND PLOTS 
%=======================================
disp("Starting Postporcess")

% Save data
save('beam_dataSolv_torsion.mat','u_hat','Fx','Fy','Fz','Mx','My','Mz','uy_modes','uz_modes','thetax_modes','freq');

% Include plot functions

% Extract the values of xn(:,1) and the selected elements of u_hat
xn_values = xn(:, 1);
u_z = u_hat(3:6:end);
thetax = u_hat(4:6:end);

% Plotting parameters
fs = 20;

% Specify the folder name
folderName = 'modes_plots_torsion';

% Create the folder if it doesn't exist
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

% Plotting the values
figure;
plot(xn_values, u_z*1e3, '-', 'LineWidth', 2, Color='red');

% Adding labels and title
xlabel('Position [m]', 'Interpreter','latex',FontSize=fs);
ylabel('Displacement $u_{z}$ [mm]','Interpreter','latex',FontSize=fs);
title('Position VS Displacement', 'Interpreter','latex',FontSize=fs);
set(gca,"ticklabelinterpreter",'latex')

% Display the grid
grid on;

fileName = fullfile(folderName, 'Displacement_uz.png');
saveas(gcf, fileName);

% Plotting the values
figure;
plot(xn_values, thetax, '-', 'LineWidth', 2, Color='red');

% Adding labels and title
xlabel('Position [m]', 'Interpreter','latex',FontSize=fs);
ylabel('Displacement $\theta_{x}$ [Rads]','Interpreter','latex',FontSize=fs);
title('Position VS Displacement', 'Interpreter','latex',FontSize=fs);
set(gca,"ticklabelinterpreter",'latex')

% Display the grid
grid on;

fileName = fullfile(folderName, 'Displacement_thetax.png');
saveas(gcf, fileName);

for i = 1:size(uz_modes,2)
    % Create a new figure for each iteration
    figure;

    % Plotting the values
    plot(xn_values, uz_modes(:,i) * 1e3, '-', 'LineWidth', 2, 'Color', 'blue');

    % Adding labels and title
    xlabel('Position [m]', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('Displacement $u_{z}$ [mm]', 'Interpreter', 'latex', 'FontSize', fs);
    title(['Mode ' num2str(i)], 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'ticklabelinterpreter', 'latex');
    
    % Display the grid
    grid on;
    
    % Save the figure in the specified folder
    fileName = fullfile(folderName, ['uz_Mode_' num2str(i) '_Plot.png']);
    saveas(gcf, fileName);
    
end

for i = 1:size(uy_modes,2)
    % Create a new figure for each iteration
    figure;

    % Plotting the values
    plot(xn_values, uy_modes(:,i) * 1e3, '-', 'LineWidth', 2, 'Color', 'green');

    % Adding labels and title
    xlabel('Position [m]', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('Displacement $u_{y}$ [mm]', 'Interpreter', 'latex', 'FontSize', fs);
    title(['Mode ' num2str(i)], 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'ticklabelinterpreter', 'latex');
    
    % Display the grid
    grid on;
    
    % Save the figure in the specified folder
    fileName = fullfile(folderName, ['uy_Mode_' num2str(i) '_Plot.png']);
    saveas(gcf, fileName);
    
end

for i = 1:size(thetax_modes,2)
    % Create a new figure for each iteration
    figure;

    % Plotting the values
    plot(xn_values, thetax_modes(:,i), '-', 'LineWidth', 2, 'Color', 'magenta');

    % Adding labels and title
    xlabel('Position [m]', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('Displacement $\theta_{x}$ [Rads]', 'Interpreter', 'latex', 'FontSize', fs);
    title(['Mode ' num2str(i)], 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'ticklabelinterpreter', 'latex');
    
    % Display the grid
    grid on;
    
    % Save the figure in the specified folder
    fileName = fullfile(folderName, ['thetax_Mode_' num2str(i) '_Plot.png']);
    saveas(gcf, fileName);
    
end

disp("Postprocess Finished Succesfully!")