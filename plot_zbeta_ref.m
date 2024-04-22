function plot_zbeta_ref()
    %%plot Z and beta2 for the ensemble mean
    
    load('reference.mat')

    msize = 150;

    scatter(reference.Z(1), reference.beta2(1), msize,'o', 'filled','c');
    scatter(reference.Z(2), reference.beta2(2), msize,'o', 'filled', 'r');
    scatter(reference.Z(3), reference.beta2(3), msize,'o', 'filled', 'm');
    errorbar(reference.Z(4), reference.beta2(4), 50,50,3.7e5,3.7e5,'o', 'Color', 'g', 'MarkerSize', 10, 'LineWidth',1.5);
    scatter(reference.Z(5), reference.beta2(5), msize,'o', 'filled', 'b');
    scatter(reference.Z(6), reference.beta2(6), msize,'o', 'filled', 'y');% 0,0,reference.Z(6)-8.75e6, 11e6-reference.Z(6), 'k');
    errorbar(reference.Z(7), reference.beta2(7), 100,100,3.9e5,3.9e5,'o', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 10, 'LineWidth',1.5);

    %scatter(Z_ensemble_mean, beta2_ensemble_mean, 'o', 'k');
    %errorbar(Z_ensemble_mean, beta2_ensemble_mean, beta2_ensemble_std, beta2_ensemble_std, Z_ensemble_std, Z_ensemble_std, 'o', 'Color',[0.5 0.5 0.5], 'LineWidth',1);

    xlabel('Acoustic impedance (kg m^{-2} s^{-1})');
    ylabel('Vs (m s^{-1})');
    
    %standard_figure;
        thesis_figure;