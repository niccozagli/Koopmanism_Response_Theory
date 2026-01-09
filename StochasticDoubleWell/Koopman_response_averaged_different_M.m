%% ====== Initialise parameters ======
% No need to run this section if stochasticdoublewell_main has already been run.
% This section is retained, so script may be run without any pre-loaded variables.

% Dynamics: double well potential
V = @(x,y) (x.^2 - 1).^2 + y.^2;
dV = @(x) [4.*(x(1).^3) - 4.*x(1), 2.*x(2)];
% We use these separate defs to plot vector field
dVx = @(x) 4.*x.*(x.^2 - 1);
dVy = @(y) 2.*y;

% Choose noise strength
sigma = 0.7;

tmax = 5*10^5;
dt_int = 0.001;
transient = 100;
tmax_tot = tmax + transient;
num_tsteps_int = round(tmax_tot/dt_int);

% ===EDMD Parameters===
% highest power in monomial dictionary 
Kmax = 10;
% dictionary size
kn = nchoosek(Kmax+2,Kmax);
% Choose EDMD timestep (flight time)
dt = 0.05;

obs_x = @(x,y) x;
obs_rad = @(x,y) x.^2 + y.^2;

% First perturbation
gamma_an1 = @(x,y) (2/sigma.^2)*(4*x.^3 - 4*x);
% Second perturbation
gamma_an2 = @(x,y) (32/sigma.^2)*(x.^4 - x.^2) - 4;

% Decompose polynomial observables by hand
h_x = zeros(kn,1); h_rad = zeros(kn,1);

h_x(3) = 1;
h_rad(4) = 1; h_rad(6) = 1;

% decompose these functions wrt monomial dictionary
h_g1 = zeros(kn,1); 
h_g2 = zeros(kn,1); 

h_g1(3) = -8/sigma.^2;
h_g1(10) = 8/sigma.^2;
h_g2(1) = -4;
h_g2(6) = -32/sigma.^2;
h_g2(15) = 32/sigma.^2;

%% ======== MAIN ========
%% Loop over routines using trajectories of different lengths as EDMD input

% Set numtrajs to 60 to reproduce results in paper
numtrajs = 10;
numlengths = 4;

% These arrays just store the number of timesteps/tmax used in each run,
% used only for truncation.
% Parameters used in integration set previously
M_full = tmax/dt;
M_arr = M_full./(10.^(0:numlengths-1));
t_arr = tmax./(10.^(0:numlengths-1));

% Now create arrays of the points at which to evaluate the GFs
tmax_odd = 50;
tmax_even = 1.2;
% resolution required:
res_odd = 1;
res_even = 0.1;

tplot_odd = 0:dt*res_odd:tmax_odd;
tplot_even = 0:dt*res_even:tmax_even;

% Finally, arrays to store all of the GFs generated
greensfns_odd = zeros(numtrajs,numlengths,length(tplot_odd));
greensfns_even = zeros(numtrajs,numlengths,length(tplot_even));

for i = 1:numtrajs
    fprintf('trajectory %d of %d',i,numtrajs)
    tic
    x0 = -1.5 + 3*rand(1,2)
    y_full = integrate2D(x0,num_tsteps_int,dV,sigma,dt_int);
    y_full = y_full(round(transient/dt_int)+1:end,:);
    filter = round(dt/dt_int);
    y = y_full(1:filter:end,:);
    toc
    lambdas = cell(1,numlengths);
    Xis = cell(1,numlengths);
    Ws = cell(1,numlengths);
    Gs = cell(1,numlengths);
    
    tic
    for j = 1:numlengths
        [Xi_temp,W_temp,lambda_temp,G_temp] = full_hermitian_EDMD_routine(y(1:M_arr(j)+1,:),Kmax);
        lambdas{j} = lambda_temp;
        Xis{j} = Xi_temp;
        Ws{j} = W_temp;
        Gs{j} = G_temp;
    
        v_g1_temp = Ws{j}'*h_g1; v_g2_temp = Ws{j}'*h_g2;

        v_x_temp = Ws{j}'*h_x; v_rad_temp = Ws{j}'*h_rad;
    
        scalarprod_matrix_temp = Xis{j}'*Gs{j}*Xis{j};

        GreenFn_x_temp = KoopmanCorrs(v_x_temp,v_g1_temp,lambdas{j},scalarprod_matrix_temp);
        GreenFn_rad_temp = KoopmanCorrs(v_rad_temp,v_g2_temp,lambdas{j},scalarprod_matrix_temp);

        greensfns_odd(i,j,:) = GreenFn_x_temp(tplot_odd/dt);
        greensfns_even(i,j,:) = GreenFn_rad_temp(tplot_even/dt);
    
    end
    toc
end

%% If more trajectories are required, re-run above, and join to form _tot arrays

greensfns_even_tot = greensfns_even;
greensfns_odd_tot = greensfns_odd;

