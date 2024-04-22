function plot_zpois_ref()
    %%plot Z and poisson for the ensemble mean
    
    load('reference.mat')

    msize = 150;


    scatter(reference.Z(1), reference.poisson(1), msize,'o', 'filled','c');
    scatter(reference.Z(2), reference.poisson(2), msize,'o', 'filled', 'r');
    scatter(reference.Z(3), reference.poisson(3), msize,'o', 'filled', 'm');
    errorbar(reference.Z(4), reference.poisson(4), 0.004,0.004,3.7e5,3.7e5,'o', 'Color', 'g', 'MarkerSize', 10, 'LineWidth',1.5);
    scatter(reference.Z(5), reference.poisson(5), msize,'o', 'filled', 'b');
    scatter(reference.Z(6), reference.poisson(6), msize,'o', 'filled', 'y');% 0,0,reference.Z(6)-8.75e6, 11e6-reference.Z(6), 'k');
    errorbar(reference.Z(7), reference.poisson(7), 0.08,0.08,3.9e5,3.9e5,'o', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 10, 'LineWidth',1.5);


    %scatter(Z_ensemble_mean, poisson_ensemble_mean, 'o', 'k');
    %errorbar(Z_ensemble_mean, poisson_ensemble_mean, poisson_ensemble_std, poisson_ensemble_std, Z_ensemble_std, Z_ensemble_std, 'o', 'Color',[0.5 0.5 0.5], 'LineWidth',1);

    xlabel('Acoustic impedance (kg m^{-2} s^{-1})');
    ylabel('Poisson ratio');
    
    standard_figure;
    ylim([0 0.5])
    pbaspect([1 1 1])