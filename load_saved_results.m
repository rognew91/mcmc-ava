addpath('~/Documents/MATLAB/crewes/syntraces');
addpath('~/Documents/MATLAB/crameri_v1.07/crameri');

BED = inputs.BED;
includePS = inputs.includePS;
includeSS = inputs.includeSS;
number_iterations = inputs.number_iterations;
burn_in = inputs.burn_in;
perturbation_scaling = inputs.perturbation_scaling;
n_ensemble = inputs.n_ensemble;
figson = inputs.figson;
synthetic = inputs.synthetic;
save_data = inputs.save_data;
anginc_PP = inputs.anginc_PP;
anginc_PS = inputs.anginc_PS;
anginc_SS = inputs.anginc_SS;
Rpp_input = inputs.Rpp_input;
Rps_input = inputs.Rps_input;
Rss_input = inputs.Rss_input;
accept_cred_interval = inputs.accept_cred_interval;
cutoffAngle = inputs.cutoffAngle;
cutoffAngle_PP = inputs.cutoffAngle_PP;
cutoffAngle_PS = inputs.cutoffAngle_PS;
cutoffAngle_SS = inputs.cutoffAngle_SS;
use_all_data = inputs.use_all_data;
sigs_PP = inputs.sigs_PP;
sigs_PS = inputs.sigs_PS;
sigs_SS = inputs.sigs_SS;

std_rho1 = priors.std_rho1;
std_alpha1 = priors.std_alpha1;
std_beta1 = priors.std_beta1;
lims_rho2 = priors.lims_rho2;
lims_alpha2 = priors.lims_alpha2;
lims_beta2 = priors.lims_beta2;
lims_poisson = priors.lims_poisson;

models_saved = results.models_saved;
log_likelihood_saved = results.log_likelihood_saved;
log_posterior_saved = results.log_posterior_saved;
log_prior_saved = results.log_prior_saved;
m_best = results.m_best;
log_posterior_best = results.log_posterior_best;
besti = results.besti;
cred_interval = results.cred_interval;
Z_best = results.Z_best;
poisson_best = results.poisson_best;
Z_saved = results.Z_saved;
poisson_saved = results.poisson_saved;
median_model = results.median_model;
iqr_model = results.iqr_model;
mode_model = results.mode_model;
models_pbi = results.models_pbi;


    x = models_saved(burn_in+1:number_iterations, 4);
    y = models_saved(burn_in+1:number_iterations, 6);
    z = models_saved(burn_in+1:number_iterations, 2);
    %c = log_likelihood_saved(burn_in+1:number_iterations);
    c = exp(log_posterior_saved(burn_in+1:number_iterations));


    %foo(1:burn_in) = 3000;
    %{
    if figson
        figure; hold on;
        %plot3(models_saved(1:burn_in, 4), models_saved(1:burn_in, 6), foo, 'r')
        plot3(x,y,z);
        cla
        patch([x; nan], [y; nan], [z; nan], [c; nan], 'EdgeColor', 'interp', 'FaceColor', 'none');
        xlabel('vp'); ylabel('vs'); zlabel('density')
        %xlim([0 5000]); ylim([0 3000]); zlim([0 3000]);
    end
    %}

    %[log_likelihood_best, besti] = max(log_likelihood_saved);
    [log_posterior_best, besti] = max(log_posterior_saved);
    m_best = models_saved(besti, :);

    %models_pbi = post burn in
    models_pbi = models_saved(burn_in+1:number_iterations, :);
    log_posterior_pbi = log_posterior_saved(burn_in+1:number_iterations);

    %find uncertainties using credible interval defined in preamble
    for err_thresh = 0.1:-0.0001:0.0001
        log_level = log(err_thresh) + log_posterior_best;
        modelnum_above_level = find(log_posterior_pbi > log_level);
        models_above_level = models_pbi(modelnum_above_level, :);
        max_of_models = max(models_above_level(:,:));
        min_of_models = min(models_above_level(:,:));
        pm_of_models = (max_of_models - min_of_models)/2;



        cred_interval = length(modelnum_above_level)/(number_iterations-burn_in);

        if cred_interval > accept_cred_interval
            break;
        end

    end

    models_pbi_sorted = sortrows(models_pbi, [4 6 2 3 5 1]);
    bar = find(models_pbi_sorted(:,4) == min_of_models(4));
    conf_isovalue = max(log_posterior_pbi(bar));



%model AVO curves of best models.
    coefPP_best = zoeppritz(m_best(1), m_best(3), m_best(5), m_best(2), m_best(4), m_best(6), 1, 1, 0, anginc_PP);
    Rpp_best = real(coefPP_best);

    coefPS_best = zoeppritz(m_best(1), m_best(3), m_best(5), m_best(2), m_best(4), m_best(6), 1, 2, 0, anginc_PS);
    Rps_best = real(coefPS_best);

    coefSS_best = zoeppritz(m_best(1), m_best(3), m_best(5), m_best(2), m_best(4), m_best(6), 2, 2, 0, anginc_SS);
    Rss_best = real(coefSS_best);

    %%find acoustic impedance and poisson's ratio
    alphabeta = m_best(4)/m_best(6);
    Z_best = m_best(4)*m_best(2);
    poisson_best = (alphabeta^2 - 2)/(2*(alphabeta^2 - 1));

    Z_saved = models_saved(burn_in+1:number_iterations,4).*models_saved(burn_in+1:number_iterations,2);
    alphabeta_saved = models_saved(burn_in+1:number_iterations,4)./models_saved(burn_in+1:number_iterations,6);
    poisson_saved = (alphabeta_saved.^2 - 2)./(2.*(alphabeta_saved.^2 - 1));

    %% find median model

    median_model = median(models_saved(burn_in+1:number_iterations,:));
    iqr_model = iqr(models_saved(burn_in+1:number_iterations,:));

    coefPP_med = zoeppritz(median_model(1), median_model(3), median_model(5), median_model(2), median_model(4), median_model(6), 1, 1, 0, anginc_PP);
    Rpp_med = real(coefPP_med);

    coefPS_med = zoeppritz(median_model(1), median_model(3), median_model(5), median_model(2), median_model(4), median_model(6), 1, 2, 0, anginc_PS);
    Rps_med = real(coefPS_med);

    coefSS_med = zoeppritz(median_model(1), median_model(3), median_model(5), median_model(2), median_model(4), median_model(6), 2, 2, 0, anginc_SS);
    Rss_med = real(coefSS_med);


    %% plost all the results
