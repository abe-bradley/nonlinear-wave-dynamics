%% Initialize notebook
clear
clc

eta0 = 376.73; % vacuum impedance in [Ohms]

%% Grid
% Time
T = 50; % in [ns]
Nt = 2^13; % keep it an even number

dt = T/Nt;
t  = (0:Nt-1)*dt;
T  = Nt*dt;

% Space
Z = 1e5; % in [um]
Nz = 2^14;

dz = Z/(Nz - 1);
z = 0:dz:Z;

% Fourier space
df = 1/T; % in [GHz]
f = df*[0:Nt/2-1, -Nt/2:-1].';
omega = 2*pi*f;

%% Input parameters in natural units
% Excitation properties
lambda0 = 1.064; % in [um]
k0 = 2*pi/lambda0;
P0 = 200; % in [kW]

% Material properties
[beta0,beta1,beta2] = sellmeier_silica(lambda0); % in [1/um], [ns/um], and [ns^2/um]
alpha = 4.6e-11; % in [1/um]

% Nonlinear properties
n2 = 2.7e-5; % in [um^2/kW]
A_eff = 100; % in [um^2]
gamma = k0*n2/A_eff; % in [um/kW]

%% Normalizations
t0 = sqrt(abs(beta2)/2/gamma/P0);
z0 = 1/gamma/P0;

s = sign(beta2);
a = alpha/2/gamma/P0;
b = sqrt(2*beta1^2/abs(beta2)/gamma/P0);

%% Initialize the field parameters
N_mem = 20;
dF = 0.1; % in [GHz]
sigma = sqrt(2*pi*dF*dt);

%% Beam propagation loop
dZ = dz/z0;
dT = dt/t0;
Omega = omega*t0;

HL = exp(dZ/2 * (-a + 1i*(s*Omega.^2 - b*Omega))); % half-step linear propagator

% Time domain
P = zeros(1, Nz);

U_in  = zeros(Nt,N_mem);
U_out = zeros(Nt,N_mem);

tic
for ii = 1:N_mem
    progress_points(ii,N_mem)

    phi = cumsum(sigma*randn(Nt,1));
    U0  = exp(1i*phi);

    P(1) = P(1) + mean(abs(U0).^2);

    U_in(:,ii) = U0;
    for mm = 1:Nz-1
        U1 = ifft(fft(U0).*HL);
        U2 = exp(1i*dZ*abs(U1).^2).*U1;
        U0 = ifft(fft(U2).*HL);

        P(mm + 1) = P(mm + 1) + mean(abs(U0).^2);
    end
    U_out(:,ii) = U0;
end
toc

% Frequency domain
Uw_in  = fft(U_in,[],1);
Uw_out = fft(U_out,[],1);

Iw_in  = mean(abs(Uw_in).^2,2);
Iw_out = mean(abs(Uw_out).^2,2);

%% Plot temporal and spectral info
% Temporal
figure
tiledlayout(2,1)
nexttile
plot(t,abs(U_in(:,1)).^2)
xlabel("Time: t [ns]")
ylabel("Initial intensity")
nexttile
plot(t,angle(U_in(:,1)),'.')
yline(0)
xlabel("Time: t [ns]")
ylabel("Initial phase")
    
% Frequency domain
figure
tiledlayout(2,1)
nexttile
plot(f,Iw_in)
xlabel("Frequency: f [GHz]")
ylabel("Initial power spectrum")
nexttile
plot(f,Iw_out)
xlabel("Frequency: f [GHz]")
ylabel("Final power spectrum")

%% Dimensionalize power then plot
P = P0*P/N_mem; % power in [kW]

figure
plot(z,P)
xlabel("Longitudinal axis: z [\mum]")
ylabel("Power: P(z) [kW]")