
%% save the results of the inversion

inputs.BED = BED;
inputs.includePS = includePS;
inputs.includeSS = includeSS;
inputs.number_iterations = number_iterations;
inputs.burn_in = burn_in;
inputs.perturbation_scaling = perturbation_scaling;
inputs.n_ensemble = n_ensemble;
inputs.figson = figson;
inputs.synthetic = synthetic;
inputs.save_data = save_data;
inputs.anginc_PP = anginc_PP;
inputs.anginc_PS = anginc_PS;
inputs.anginc_SS = anginc_SS;
inputs.Rpp_input = Rpp_input;
inputs.Rps_input = Rps_input;
inputs.Rss_input = Rss_input;
inputs.accept_cred_interval = accept_cred_interval;
inputs.cutoffAngle = cutoffAngle;
inputs.cutoffAngle_PP = cutoffAngle_PP;
inputs.cutoffAngle_PS = cutoffAngle_PS;
inputs.cutoffAngle_SS = cutoffAngle_SS;
inputs.use_all_data = use_all_data;
inputs.sigs_PP = sigs_PP;
inputs.sigs_PS = sigs_PS;
inputs.sigs_SS = sigs_SS;
inputs.PPdatafile = PPdatafile;
inputs.PSdatafile = PSdatafile;

priors.std_rho1 = std_rho1;
priors.std_alpha1 = std_alpha1;
priors.std_beta1 = std_beta1;
priors.lims_rho2 = lims_rho2;
priors.lims_alpha2 = lims_alpha2;
priors.lims_beta2 = lims_beta2;
priors.lims_poisson = lims_poisson;

results.models_saved = models_saved;
results.log_likelihood_saved = log_likelihood_saved;
results.log_posterior_saved = log_posterior_saved;
results.log_prior_saved = log_prior_saved;
results.m_best = m_best;
results.log_posterior_best = log_posterior_best;
results.besti = besti;
results.cred_interval = cred_interval;
results.Z_best = Z_best;
results.poisson_best = poisson_best;
results.Z_saved = Z_saved;
results.poisson_saved = poisson_saved;
results.median_model = median_model;
results.iqr_model = iqr_model;
results.mode_model = mode_model;
results.models_pbi = models_pbi;
results.elapsed_time = toc;
results.misfit_saved = misfit_saved;
results.runningmedianbeta = runningmedianbeta;
results.runningmedianZ = runningmedianZ;
results.runningmedianpois = runningmedianpois;

%% get time and folder name for saving
time = datestr(now, 'yyyy-mm-dd-hh-MM');
formatOut = 'yyyy-mm-dd';

if ~includePS
    label1 = 'P';
elseif ~includeSS && includePP
        label1='C';
elseif ~includeSS && ~includePP
        label1 = 'c';
else 
    label1='S';
end

if use_all_data
    label2 = 'AL';
else
    label2 = num2str(cutoffAngle_PP);
end
label = strcat(label1,label2);

invfoldername = strcat(foldername,'/',label,'_',time);
figfoldername = strcat(invfoldername, '/figures');

savedatafile = strcat(invfoldername,'/',label,'_',time,'.mat');
if ~exist(foldername, 'dir')
    mkdir(foldername)
end
if ~exist(invfoldername, 'dir')
    mkdir(invfoldername)
end

save(savedatafile, 'inputs', 'results', 'priors');
fprintf('\nsaving into file\n %s\n\n', savedatafile);

if save_figs
    if ~exist(figfoldername, 'dir')
    mkdir(figfoldername)
end

    FigList = findobj(allchild(0), 'flat','Type','figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
      savefig(FigHandle, fullfile(figfoldername, strcat(FigName, '.fig')));
    end
    fprintf('\nsaving figs into folder %s\n\n', figfoldername);

end
