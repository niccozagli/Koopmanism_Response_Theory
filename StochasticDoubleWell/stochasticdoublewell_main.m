% ====== INITIALIZATION ======

% Dynamics: double well potential
V = @(x,y) (x.^2 - 1).^2 + y.^2;
% & its gradient
dV = @(x) [4.*(x(1).^3) - 4.*x(1), 2.*x(2)];
% We use these separate defs to plot vector field
dVx = @(x) 4.*x.*(x.^2 - 1);
dVy = @(y) 2.*y;

% Choose noise strength
sigma = 0.7;

% Choose initial cond. (x0,y0) randomly within square 
x0 = -1.5 + 3*rand(1,2);

% integration parameters
tmax = 5*10^5;
dt_int = 0.001;
transient = 100;
tmax_tot = tmax + transient;
num_tsteps_int = round(tmax_tot/dt_int);

% ===EDMD Parameters===
% Choose highest power in monomial dictionary 
Kmax = 10;
% dictionary size
kn = nchoosek(Kmax+2,Kmax);
% Choose EDMD timestep (flight time)
dt = 0.05;

%% ===Discrete grids for plotting, numerically decomposing functions ===
% Define granularity
xdiscrete = -2:0.05:2; ydiscrete = xdiscrete;
% Create grid
[xgrid,ygrid] = meshgrid(xdiscrete,ydiscrete);

% Count bins
nbins_x = length(xdiscrete)-1;nbins_y = length(ydiscrete)-1;
nbins_tot = nbins_x*nbins_y;

% Create grid at midpoints for plotting mu
xdisc_mid = zeros(nbins_x,1);
ydisc_mid = zeros(nbins_y,1);
for i=1:nbins_x
    xdisc_mid(i) = (xdiscrete(i+1) + xdiscrete(i))/2;
end
for i=1:nbins_y
    ydisc_mid(i) = (ydiscrete(i+1) + ydiscrete(i))/2;
end
[xgrid_mid,ygrid_mid] = meshgrid(xdisc_mid,ydisc_mid);

%% =============== MAIN ===============

% Integrate
y_full = integrate2D(x0,num_tsteps_int,dV,sigma,dt_int);

% first remove transient
y_full = y_full(round(transient/dt_int)+1:end,:);

% then filter trajectory to use for EDMD
filter = round(dt/dt_int);
y = y_full(1:filter:end,:);

disp('Integration completed')

% === Get measure probabilistically ===
mu = histcounts2(y_full(:,1),y_full(:,2),xdiscrete,ydiscrete,'Normalization','pdf');
mu = mu.';
disp('Invariant density updated')

% ====== EDMD ======
% Monomial dictionary is hard-coded into function
[G,A] = EDMD_2D_vectorised(y,Kmax);

% Hermitian DMD: comment these lines to turn off:
A_h = (A+A')/2;
A = A_h;
% This produces virtually identical results to EDMD
 
K = pinv(G)*A;

% Get spectral properties of K
[Xi,W,lambda] = get_spectral_properties(K);
disp('EDMD completed and spectral properties calculated and updated')

% eigenvalues of generator: decay rates
rates = log(lambda)/dt;

%% ====== ANALYSIS ======

% Reconstruct observables:
obs_x = @(x,y) x;
obs_y = @(x,y) y;
obs_rad = @(x,y) x.^2 + y.^2;
obs_3 = @(x,y) cos(2*x) + sin(2*y);
obs_4 = @(x,y) sin(2*x) + cos(2*y);

% Dictionary coefficients of polynomial observables by hand
h_x = zeros(kn,1); h_y = zeros(kn,1); h_rad = zeros(kn,1);
h_y(2) = 1;
h_x(3) = 1;
h_rad(4) = 1; h_rad(6) = 1;

% Get Koopman coefficients
v_y = W'*h_y; v_x = W'*h_x; v_rad = W'*h_rad;

% Decompose non-polynomial observables
% Symbolic integration is also possible, but slower
obsnum_3 = obs_3(xgrid,ygrid);
[hnum_3, vnum_3] = decompose_genfn_numerical(xdiscrete,ydiscrete,obsnum_3,Kmax,W);

obsnum_4 = obs_4(xgrid,ygrid);
[hnum_4, vnum_4] = decompose_genfn_numerical(xdiscrete,ydiscrete,obsnum_4,Kmax,W);

% === Compute correlation functions of observables ===
% Koopman correlations
scalarprod_matrix= Xi'*G*Xi;
% Odd functions
Kcorrfn_x = KoopmanCorrs(v_x,v_x,lambda,scalarprod_matrix);
Kcorrfn_4 = KoopmanCorrs(vnum_4,vnum_4,lambda,scalarprod_matrix);
% Even functions
Kcorrfn_y = KoopmanCorrs(v_y,v_y,lambda,scalarprod_matrix);
Kcorrfn_rad = KoopmanCorrs(v_rad,v_rad,lambda,scalarprod_matrix);
Kcorrfn_3 = KoopmanCorrs(vnum_3,vnum_3,lambda,scalarprod_matrix);

% ====== Response functions ======
% Define Gamma analytically & decompose
% First perturbation
gamma_an1 = @(x,y) (2/sigma.^2)*(4*x.^3 - 4*x);
% Second perturbation
gamma_an2 = @(x,y) (32/sigma.^2)*(x.^4 - x.^2) - 4;

% decompose these functions wrt monomial dictionary
h_g1 = zeros(kn,1); 
h_g2 = zeros(kn,1); 

h_g1(3) = -8/sigma.^2;
h_g1(10) = 8/sigma.^2;
h_g2(1) = -4;
h_g2(6) = -32/sigma.^2;
h_g2(15) = 32/sigma.^2;

% Koopman modes
v_g1 = W'*h_g1; v_g2 = W'*h_g2;

% Calculate response functions for prev. defined observables to these perturbations
GreenFn_x = KoopmanCorrs(v_x,v_g1,lambda,scalarprod_matrix);
GreenFn_rad = KoopmanCorrs(v_rad,v_g2,lambda,scalarprod_matrix);

%% Compute correlation functions numerically on independent trajectory

% First generate new, independent trajectory 
tmax_tot2 = 5*10^5;
transient2 = 10^3;
tmax_2 = tmax_tot2 + transient2;
num_tsteps_int_2 = round(tmax_2/dt_int);
x0_2 = -1.5 + 3*rand(1,2);

% Overwrite original trjactory to save memory
y_full = integrate2D(x0_2,num_tsteps_int_2,dV,sigma,dt_int);

dt_2 = 0.01;
filter_2 = round(dt_2/dt_int);

y_2 = y_full(round(transient2/dt_int)+1:filter_2:end,:);

% Numerical correlations
max_lag = 200/dt_2;
% Odd functions (in x)
% Only need to evaluate lags once
[cf_x,lags] = get_cf_2D(y_2,obs_x,obs_x,max_lag);
cf_4 = get_cf_2D(y_2,obs_4,obs_4,max_lag);
% Even functions
cf_y = get_cf_2D(y_2,obs_y,obs_y,max_lag);
cf_rad = get_cf_2D(y_2,obs_rad,obs_rad,max_lag);
cf_3 = get_cf_2D(y_2,obs_3,obs_3,max_lag);

