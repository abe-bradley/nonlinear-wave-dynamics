%% Initialize notebook
clear
clc

eta0 = 376.73; % vacuum impedance in [Ohms]

%% Input parameters in natural units
% Choice of domain size and resolution
lx = 40; % in [um]
lz = 400; % in [um]

Nx = 10^3 + 1;
Nz = 10^4 + 1;

% Material properties
n0 = 1.5; % reference, for normalization
n_co = 3.332; % core
n_cl = 1.444; % cladding

n2 = 2e-13; % in [cm^2/W]
n2_bar = n0/2/eta0*n2 * 1e8; % in [um^2/V^2]

% Vacuum wavelength
lambda = 1.55; % in [um]
k0 = 2*pi/lambda; % in [1/um]

%% Scaling for spatial dimensions
x0 = 1; % our choice, in [um]
z0 = 2*k0*n0*x0^2; % in [um]

%% Dimensionless spatial grids
% Transverse direction
Lx = lx/x0;
DX = Lx/(Nx - 1);
X = 0:DX:Lx;
Xs = X - Lx/2; % centered-coordinate shift

% Paraxial longitudinal direction
Lz = lz/z0;
DZ = Lz/(Nz - 1);
Z = 0:DZ:Lz;

% Fourier space of the transverse direction
Kx = 2*pi/Lx*[0:ceil(Nx/2)-1 -floor(Nx/2):-1]';

%% Dimensional spatial grids
z = Z*z0;
x = X*x0;
dz = DZ*z0;
dx = DX*x0;
xs = Xs*x0;
kx = Kx/x0; % spatial frequency in [1/um]

%% Initialize index variation and field
U = zeros(Nx,Nz); % normalized field
nL = n_cl*ones(Nx,Nz); % cladding refractive index profile

%% Generate structure
A = 4; % dimensionless waveguide radius
Na = floor(A/DX);
wg = floor(Nx/2) + (-Na:Na);
nL(wg,:) = n_co;

%% Example structures and excitations (optional)
% To use an example, comment out the default Generate structure block and the
% Automatic (mode) excitation block, then uncomment one example below.
% All lengths below are dimensionless; C sets the normalized amplitude.

% % One soliton
% nL(:) = n_cl;
% C = 1; % normalized amplitude (integer order of soliton)
% W = sqrt(2);
% U(:,1) = C*sech(Xs/W);
% % U(:,1) = C*exp(-(Xs/W/1.7).^2);

% % Two soliton
% nL(:) = n_cl;
% B = 4; % half distance between solitons
% phi = pi/2; % second soliton phase shift
% C = 1; % normalized amplitude (integer order of soliton)
% W = sqrt(2);
% U(:,1) = C*(sech((Xs - B)/W) + exp(1i*phi)*sech((Xs + B)/W));

% % Waveguide
% nL(:) = n_cl;
% A = 3; % waveguide radius
% Na = floor(A/DX);
% wg = floor(Nx/2) + (-Na:Na);
% nL(wg,:) = n_co;
%
% C = 1e-4; % normalized amplitude
% W = sqrt(13);
% U(:,1) = C*exp(-(Xs/W).^2);

% % Coupler
% nL(:) = n_cl;
% A = 3; % waveguide radii
% B = 0.1; % distance between waveguides
% Na = floor(2*A/DX);
% Nb = floor(B/DX);
% wg1 = floor(Nx/2) + Nb + (0:Na);
% wg2 = floor(Nx/2) - Nb - (0:Na);
% nL([wg1,wg2],:) = n_co;
%
% C = 5e-1; % normalized amplitude
% W = sqrt(13);
% U(:,1) = C*exp(-((Xs - (A + B))/W).^2);

%% Mode solver
% Solve the eigenproblem for the Helmholtz operator in the initial plane
tic
H = laplacian_1d(Nx) + spdiags((k0*dx*nL(:,1)).^2, 0, Nx, Nx);
[Q,D] = eig(H);
toc

% Physical propagation constants and modes normalized on the X grid
beta = sqrt(diag(D))/dx; % in [1/um]
n_eff = beta/k0;
Psi = Q/sqrt(DX); % sum(abs(Psi).^2,1)*DX = 1

% Sort modes from highest effective refractive index to lowest
n_eff = flip(n_eff);
Psi = flip(Psi,2);
psi = Psi/sqrt(x0); % modes normalized on the physical x grid, in [1/sqrt(um)]

% Number of guided modes
M = sum((imag(n_eff) == 0) & (real(n_eff) > n_cl));

% Plot
figure
plot(n_eff(1:M),'o')
xlabel("Mode index")
ylabel('Effective refractive index: n_eff')

%% Manual excitation (optional)
% Alternative to the examples above; disable automatic mode excitation.
% C = 0.85; % normalized amplitude
% W = sqrt(2.5);
% U(:,1) = C*exp(-(Xs/W/sqrt(pi)).^2);
% % U(:,1) = C*sech(Xs/W);

%% Automatic (mode) excitation
J = zeros(Nx,Nz);

P0 = 1; % normalized excitation power [a.u.]
N = 2;
J(1:N,1) = sqrt(P0/N);

U(:,1) = Psi*J(:,1);

%% Normalize nonlinear and potential terms
delta = 1/(k0*x0)^2; % potential
u0 = 1./sqrt(2*nL*n2_bar)/k0/x0; % field scale in [V/um], note its spatial dependence

V = (nL.^2 - n0^2)/delta;

%% Plot the given conditions
% Initial field profile
figure
plot(x,abs(U(:,1)).^2)
xlabel("Transverse axis: x [\mum]")
ylabel('Initial excitation: U(z=0) [a.u.]')

% Index profile
figure
imagesc(z,x,nL)
xlabel("Longitudinal axis: z [\mum]")
ylabel("Transverse axis: x [\mum]")
c = colorbar;
c.Label.String = "Refractive index: n_L(x,z)";

%% Beam propagation loop
tic
for mm = 1:Nz-1
    U0 = U(:,mm);
    U1 = ifft(fft(U0).*exp(-1i*DZ/2*Kx.^2));       % First step (linear, half-step propagation)
    U2 = exp(1i*DZ*(V(:,mm) + abs(U1).^2)).*U1;    % Second step (nonlinear, inhomogeneous propagation)
    U3 = ifft(fft(U2).*exp(-1i*DZ/2*Kx.^2));       % Third step (linear, half-step propagation)

    U(:,mm+1) = U3;
end
toc

u = U.*u0; % dimensionalize the field

%% Field intensity
I = nL/2/eta0.*abs(u).^2; % intensity in [W/um^2]

figure
imagesc(z,x,I)
colormap(ot_colormap('jet'))
xlabel("Longitudinal axis: z [\mum]")
ylabel("Transverse axis: x [\mum]")
c = colorbar;
c.Label.String = "Intensity: I(x,z) [W/um^2]";

%% Evolution of power
% The integral over x gives power per unit length in y, in [W/um].
% To report total power, assume an effective width in y equal to x0.
% This width is an additional physical assumption about the 1+1D problem.
P = x0*sum(I*dx); % power in [W]

figure
plot(z,P)
xlabel("Longitudinal axis: z [\mum]")
ylabel("Power: P(z) [W]")
