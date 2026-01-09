% This script uses arrays of GFs generated in Koopman_response_averaged_different_M.m
% to plot average and variance of response functions genrated using EDMD with trajectories of differing lengths
% These are compared with experimentally-obtained GFs from response_experiments_doublewell.m

% First: colours used in plots
colarray = ['b', 'y', 'g','r'];

%% Loop over M, plotting average GFs

count = 1;

for i = flip(1:numlengths)

    mean_data_odd = squeeze(mean(greensfns_odd_tot(:,i,:)));
    mean_data_even = squeeze(mean(greensfns_even_tot(:,i,:)));

    std_data_odd = squeeze(std(greensfns_odd_tot(:,i,:)));
    std_data_even = squeeze(std(greensfns_even_tot(:,i,:)));

    upper_bd_even = mean_data_even + std_data_even;
    lower_bd_even = mean_data_even - std_data_even;

    upper_bd_odd = mean_data_odd + std_data_odd;
    lower_bd_odd = mean_data_odd - std_data_odd;

    figure(1)
    box on
    hold on
    % First plot \pm 1 std as shaded region
    fill([tplot_odd, fliplr(tplot_odd)], [upper_bd_odd', fliplr(lower_bd_odd')],colarray(count),'FaceAlpha', 0.05+0.1*count, 'EdgeColor', 'none','HandleVisibility','off');
    % Then plot average Koopman GF as coloured line
    plot(tplot_odd,mean_data_odd,colarray(count),'LineWidth',2,'DisplayName',sprintf('$T = %g \\times 10^{%d}$', t_arr(i)/10^floor(log10(t_arr(i))), floor(log10(t_arr(i)))))

    ylim([-0.05 1.05])
    xlabel('$t$','Interpreter','latex')
    ylabel('$G_{f}(t)$','Interpreter','latex')
    legend('Location','northwest', 'Interpreter', 'latex')

    figure(2)
    box on
    hold on
    % First plot \pm 1 std as shaded region
    fill([tplot_even, fliplr(tplot_even)], [upper_bd_even', fliplr(lower_bd_even')],colarray(count),'FaceAlpha', 0.05+0.1*count, 'EdgeColor', 'none','HandleVisibility','off');
    % Then plot average Koopman GF as coloured line
    plot(tplot_even,mean_data_even,colarray(count),'LineWidth',2,'DisplayName',sprintf('$T = %g \\times 10^{%d}$', t_arr(i)/10^floor(log10(t_arr(i))), floor(log10(t_arr(i)))))

    xlabel('$t$','Interpreter','latex')
    ylabel('$G_{f}(t)$','Interpreter','latex')
    legend('Interpreter', 'latex')

    count = count+1;
end


%% Add numerical response experiments to plot. 
% These should be loaded by running script response_experiments_doublewell.m
figure(1)
plot(tt_resp,gf_num_x,'k.','LineWidth',2,'MarkerSize',5,'DisplayName','Numerical exp.')
figure(2)
plot(tt_resp2,gf_num_rad,'k.','LineWidth',2,'MarkerSize',5,'DisplayName','Numerical exp.')

%% Add zoomed inset to first plot

figure(1)

% Inset axes position
inset_pos = [0.5, 0.5, 0.38, 0.38];  % [x, y, width, height] in normalized units
inset = axes('Position', inset_pos);

% Plot zoomed region in the inset
hold(inset, 'on')
box on
count = 1;
for i = flip(1:numlengths)

    mean_data_odd = squeeze(mean(greensfns_odd_tot(:,i,:)));
    mean_data_even = squeeze(mean(greensfns_even_tot(:,i,:)));

    std_data_odd = squeeze(std(greensfns_odd_tot(:,i,:)));
    std_data_even = squeeze(std(greensfns_even_tot(:,i,:)));

    upper_bd_even = mean_data_even + std_data_even;
    lower_bd_even = mean_data_even - std_data_even;

    upper_bd_odd = mean_data_odd + std_data_odd;
    lower_bd_odd = mean_data_odd - std_data_odd;


    fill([tplot_odd, fliplr(tplot_odd)], [upper_bd_odd', fliplr(lower_bd_odd')],colarray(count),'FaceAlpha', 0.1+0.1*count, 'EdgeColor', 'none','HandleVisibility','off');
    plot(tplot_odd,mean_data_odd,colarray(count),'LineWidth',2,'HandleVisibility','off')

    xlabel('$t$','Interpreter','latex')
    ylabel('$G_{f}(t)$','Interpreter','latex')

    count = count+1;

end

% Again, add numerical GF
plot(tt_resp,gf_num_x,'k.','LineWidth',2,'HandleVisibility','off')
hold(inset, 'off')

xlim(inset,[0 2])
ylim(inset,[-0.05 1.05])