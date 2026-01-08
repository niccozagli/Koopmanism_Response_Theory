%% Numerical response experiments
% Long trajectory required as input: get this, along with gradient of
% potential and noise strength, from stochasticdoublewell_main.m
%% First perturbation
pert1 = @(x) [1,0];

tmax_resp = 50;
dt_resp = 0.01;
num_tsteps_resp = round(tmax_resp/dt_resp);
tt_resp = 0:dt_resp:tmax_resp;

num_samples = 5*10^6;

% Choose strength of perturbation
eps = 0.01;
tic
gf_num_x = get_response_numerical_single_eps(obs_x,pert1,eps,num_samples,y_2,num_tsteps_resp,dt_resp,dV,sigma);
toc

%% Second perturbation
pert2 = @(x) [4.*x(:,1),zeros(size(x,1),1)];

tmax_resp2 = 1.2;
dt_resp2 = 0.005;
num_tsteps_resp2 = round(tmax_resp2/dt_resp2);
tt_resp2 = 0:dt_resp2:tmax_resp2;

num_samples = 5*10^6;

% Choose strength of perturbation
eps = 0.01;
tic
gf_num_rad = get_response_numerical_single_eps(obs_rad,pert2,eps,num_samples,y_2,num_tsteps_resp2,dt_resp2,dV,sigma);
toc