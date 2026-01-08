function greenfn_single_eps = get_response_numerical_single_eps(obs,pert,eps,num_samples,inv_meas,M,dt,dV,noise)
    % M is the number of integration timesteps in the response experiments

    % Sample the data distributed wrt the inv. meas uniformly
    filter_for_inits = round(length(inv_meas(:,1))/num_samples);
    init_conds = inv_meas(1:filter_for_inits:end,:);
    
    % Apply the perturbation to these sampled points
    inits_pert = init_conds + eps.*pert(init_conds);
    
    % Evaluate the chosen observable on the unperturbed trajectory 
    obs_unpert = obs(inv_meas(:,1),inv_meas(:,2));
    obs_unpert_av = mean(obs_unpert);

    % Perform experiments to get observable on perturbed trajectory
    obs_ens_av = average_pert_obs(obs,inits_pert,M,dV,noise,dt);
    
    greenfn_single_eps = (obs_ens_av-obs_unpert_av)/eps;

end

function obs_ens_av = average_pert_obs(obs,inits_pert,M,dV,noise,dt)
    
    % Given a set of perturbed init. conds, integrate over ensemble of trajectories, 
    % and get the average of the observable evaluated on these trajectories
    
    obs_ens_av = zeros(M+1,1);

    for init = 1:length(inits_pert)
        y_temp = integrate2D(inits_pert(init,:),M,dV,noise,dt);
        obs_temp = obs(y_temp(:,1),y_temp(:,2));
        obs_ens_av = obs_ens_av + obs_temp;
    end
    
    obs_ens_av = obs_ens_av/length(inits_pert);

end