clear
clc

%% Lattice definition
folder_path = fullfile(codes_root(),'MATLAB','Data', ...
    'System configurations','Lattice configurations');
% file_name = 'Jx_4x4';
% file_name = 'triple_test';
file_name = 'Jx_4x4';
file_path = fullfile(folder_path,file_name);

[H,B] = lattice_config(true,file_path);
% [H,B] = lattice_config(true);

M = length(H);
N_sub = length(B);

sysk = eigk(false,H,B); % Kerr properties

%% Distances
zi = 0;
zf = 100;
Nz = 1000;
z = linspace(zi,zf,Nz);

N_end = 600; % number of late propagation points used in the ensemble
I_end = Nz - N_end + 1:Nz;

%% Initial condition parameters
% Ensemble parameters
N_init = 25; % number of initial conditions
SNR = 50; % noise

% Excitation parameters
P0 = {8}; % sublattice powers
exc = {1:4}; % local mode indices excited

%% Intitial conditions for subsystem eigenstates
% Occupations
J0.sub = zeros(M,1);
for ss = 1:N_sub
    J0.sub(sysk.B{ss}(exc{ss})) = P0{ss}/length(exc{ss}); % occupations
end

% Complex amplitudes
c0.sub = zeros(M,N_init);
for nn = 1:N_init
    c0.sub(:,nn) = noiser(sqrt(J0.sub),SNR,true);
end

% Basis transformation to site-local space
a0 = sysk.sub.psi*c0.sub;

% Full-system modal representation
c0.full = sysk.full.psi\a0;
J0.full = mean(abs(c0.full).^2,2);

%% Propagation loop
a = zeros(M,Nz,N_init);

opts = odeset('RelTol',1e-12,'AbsTol',1e-12); % ODE tolerance options

progress_fig = progress_bar(0,N_init,[]);
for nn = 1:N_init
    [~,am] = ode45(@(Z,a) dnlse(Z,a,H), z, a0(:,nn), opts);
    a(:,:,nn) = am.';

    progress_fig = progress_bar(nn,N_init,progress_fig);
end

%% Modal basis transformation
c.full = reshape(sysk.full.psi\reshape(a,M,[]),M,Nz,N_init);
c.sub = reshape(sysk.sub.psi\reshape(a,M,[]),M,Nz,N_init);

J.full = abs(c.full).^2;
J.sub = abs(c.sub).^2;

%% Conserved quantities
[U,P] = kerr_conservation(true,z,c,sysk);

U_ave_lin = mean(U.full.lin.exact,2);
P_ave = mean(P.full,2);

%% Plot modal dynamics
plot_modal_dynamics(N_end,z,c,sysk)

%% Compute modal statistics and plot against thermodynamic predictions
N_membs = N_init*N_end;

statistics.full = modal_statistics(J.full(:,I_end,:));
statistics.sub = repmat(statistics.full,0,1);
for ss = 1:N_sub
    statistics.sub(ss) = modal_statistics(J.sub(sysk.B{ss},I_end,:));
end

plot_thermodynamics(statistics,U,P,sysk)
