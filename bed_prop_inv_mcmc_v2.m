%bed_prop_inv_mcmc_v2.m
%performs bayesian mcmc inversion of glaciological ava data



%inputs and outputs are in .mat files
%input data must be in the format of 1xn matrices
%n is the no of points/traces in the AVA data

%you will need the statistics and machine learning toolbox
%you will also need the CREWES zoeppritz function:https://www.crewes.org/ResearchLinks/FreeSoftware/index.php

addpath('synth_inputs')%where the synthetic ava data are
addpath('C:\Users\rognew91\OneDrive - NERC\Documents\MATLAB\crewes\syntraces');%add path of CREWES zoeppritz function

clear

TRUE = 1;
FALSE = 0;


%% inversion options
%--for editing
tic

%use BED PERMA BASEM DILAT DILATN STIFF STIFFN WATER LITH REAL
BED = 'STIFFN'; %choose bed, can use various synthetics or set BED = 'REAL' for real data
includePP = TRUE; %include PP AVA data
includePS = TRUE; %include PP AVA data
includeSS = FALSE; %include PP AVA data
poisson_prior = TRUE; %use if you want to set flat priors of poissons' ratio between 0 and 0.5
forcepol = FALSE;%force a polarity reversal at cangle_forced +- cangle_std
cangle_forced = 57;
cangle_std = 4;
noise = FALSE;%if you want to add gaussian noise to input/synth data
sigs = 0.2;%error on ava data, will be overwritten if you use real data
use_all_data = FALSE; %TRUE if you want to use all avo data, FALSE if you want to truncate at a given angle
cutoffAngle = 60; %max angle of ava data to use, set = NaN to use all points read
number_iterations = 200000; % number of iterations
perturbation_scaling = 20; %scaling for mcmc perturbation, suggest 20 for synthetics
thinning = 1; %factor to thin posteriors by
%where to find the real data:
PPdatafile = 'FILENAME-FOR-REAL-PP-DATA.mat'; %input in the format anginc_PP, Rpp, Rpp_err (vectors of angle of incidence, P reflectivity and error)
PSdatafile = 'FILENAME-FOR-REAL-PS-DATA.mat'; %input in the format anginc_PS, Rps, Rps_err (vectors)
foldername = strcat('inversion_results/',BED); %where you want to save results
n_ensemble = 1;
startscale = 200;

%% saving options
figson = TRUE;
save_data = TRUE;
save_figs = FALSE;

%%
if ~strcmp(BED, 'REAL')
    synthetic = TRUE;
    if strcmp(BED, 'PERMA')
        load 'InputSynthetic_perma.mat'
    elseif strcmp(BED, 'BASEM')
        load 'InputSynthetic_basem.mat'
        cutoffAngle_SS = 20;
    elseif strcmp(BED, 'DILAT')
        load 'InputSynthetic_dilat.mat'
        cutoffAngle_SS = 28;
    elseif strcmp(BED, 'DILATN')
        load 'InputSynthetic_dilat_n.mat'
        if includeSS == TRUE
            fprintf('Rss not defined\n');
        end
        cutoffAngle_SS = 28;
    elseif strcmp(BED, 'WATER')
        load 'InputSynthetic_water.mat'
        cutoffAngle_SS = 29;
    elseif strcmp(BED, 'LITH')
        load 'InputSynthetic_lith.mat'
        cutoffAngle_SS = 28;
    elseif strcmp(BED, 'STIFF')
        load 'InputSynthetic_stiff.mat'
        cutoffAngle_SS = 27;
    elseif strcmp(BED, 'STIFFN')
        load 'InputSynthetic_stiff_n.mat'
        if includeSS == TRUE
            fprintf('Rss not defined\n');
        end
        cutoffAngle_SS = 28;
    end
    anginc_PP = anginc;
    anginc_PS = anginc;
    anginc_SS = NaN;
    %Rss_input = NaN;
    %Rss_input = Rss;
    anginc_SS = anginc;
end

if strcmp(BED, 'REAL')
    synthetic = FALSE;
    load(PPdatafile);
    Rpp_input = 2*Rpp;
    sigs_PP = Rpp_err;
    Rpp_input(128) = NaN;
    sigs_PP(128) = NaN;
    load(PSdatafile);  
    Rps_input = Rps;
    sigs_PS = Rps_err;
    Rss_input = NaN;
    anginc_SS = NaN;
        cutoffAngle_SS = NaN;

        if ~includePP
            Rpp_input(1:length(Rpp_input))=NaN;
        end

  
end
if noise
    for i=1:length(Rpp_input)        
        Rpp_input(i) = Rpp_input(i)+normrnd(0,sigs);
    end
    for i=1:length(Rps_input)
        Rps_input(i) = Rps_input(i)+normrnd(0,sigs);
    end
end

if save_figs
    close all
end

%% set the inversion up, edit this bit
%--for editing

%data = Rpp_input; %this is the data vector
% zoeppritz takes the role of g(m)

printit = number_iterations/10;
M = 6; % numel in model
models_saved = zeros(number_iterations, M);
log_likelihood_saved = zeros(number_iterations, 1);
log_posterior_saved = zeros(number_iterations, 1);
log_prior_saved = zeros(number_iterations, 1);
n_accept = 0; % number of accepted models
n_reject = 0; % number of rejected models
burn_in = 10000;

accept_cred_interval = 0.89;%0.89;


cutoffAngle_PP = cutoffAngle;
cutoffAngle_PS = cutoffAngle;

if cutoffAngle_SS > cutoffAngle
    cutoffAngle_SS = cutoffAngle;
end

if ~strcmp(BED, 'REAL')
    sigs_PP(1:length(anginc_PP)) = sigs;
    sigs_PS(1:length(anginc_PS)) = sigs;
end
sigs_SS(1:length(anginc_SS)) = sigs;


%% define starting model, can edit this bit
%------ define the starting model ------

%ice properties, please edit do not assume!!
rho1 = 920; % in Kg/m3
alpha1 = 3810; % in m/s
beta1 = 1860;

if strcmp(BED, 'REAL')
    rho1 = 920;
    alpha1 = 3830;
    beta1 = 1906;
end

rho2 = 2000 + normrnd(0, startscale);
alpha2 = 3000 + normrnd(0, startscale);
beta2 = 1500 + normrnd(0, startscale);

m_current = [rho1, rho2, alpha1, alpha2, beta1, beta2];

%% define priors, edit this bit
%------- define priors ------------------
%priors on ice properties, remember to edit!!
std_rho1 = 20;
std_alpha1 = 20;
std_beta1 = 20;

if strcmp(BED,'REAL')
    std_rho1 = 50;
    std_alpha1 = 50;
    std_beta1 = 50;
end

% set uniform priors for bed properties, min and max of uniform dist
lims_rho2 = [920, 4000];
lims_alpha2 = [0, 8000];
lims_beta2 = [0, 5000];

lims_poisson = [0 0.5];

%xscale = 1:1:5000;
%gaussian priors centred on 'ice value' for rho1, alpha1, beta1
%prior_rho1 = (1/(sqrt(2*pi)*std_rho1))*exp(-((xscale-rho1).^2)/(2.*std_rho1.^2));
%prior_alpha1 = (1/(sqrt(2*pi)*std_alpha1))*exp(-((xscale-alpha1).^2)/(2.*std_alpha1.^2));
%prior_beta1 = (1/(sqrt(2*pi)*std_beta))*exp(-((xscale-beta1).^2)/(2.*std_beta1.^2));

fprintf('\n');

%%
%%% check and edit everything above this line to set up the inversion %%%

%% restrict avo angles
if ~use_all_data
    if ~isnan(cutoffAngle)
        %cutoffIndex = find(anginc > cutoffAngle, 1, 'first') - 1;


        cutoffIndex_PP = find(anginc_PP > cutoffAngle_PP, 1, 'first') - 1;
        cutoffIndex_PS = find(anginc_PS > cutoffAngle_PS, 1, 'first') - 1;
        cutoffIndex_SS = find(anginc_SS > cutoffAngle_SS, 1, 'first') - 1;


        anginc_PP = anginc_PP(1:cutoffIndex_PP);
        anginc_PS = anginc_PS(1:cutoffIndex_PS);
        anginc_SS = anginc_SS(1:cutoffIndex_SS);

        Rpp_input = Rpp_input(1:cutoffIndex_PP);
        Rps_input = Rps_input(1:cutoffIndex_PS);
        Rss_input = Rss_input(1:cutoffIndex_SS);

        sigs_PP = sigs_PP(1:cutoffIndex_PP);
        sigs_PS = sigs_PS(1:cutoffIndex_PS);
        sigs_SS = sigs_SS(1:cutoffIndex_SS);

    end
else
    anginc_PP = anginc_PP;
    anginc_PS = anginc_PS;
    anginc_SS = anginc_SS;
end
    
%% do the inversion
for k=1:n_ensemble
    if n_ensemble > 1
        fprintf('---sim %d---\n',k)
    end


%% do an individual inversion, plot and save results
    for i = 1:number_iterations
    %%
    
     if mod(i, printit) == 0
         fprintf('\n---- inversion in progress, model %d ----', i)
    end
    % propose a model
        m_trial = m_current;
        
        for p=1:6
            perturbation(p) = normrnd(0, perturbation_scaling);
        end
        m_trial = m_trial + perturbation;
        
    % calculate Rpp for current model
        coefPP_trial = zoeppritz(m_trial(1), m_trial(3), m_trial(5), m_trial(2), m_trial(4), m_trial(6), 1, 1, 0, anginc_PP);
        Rpp_trial = real(coefPP_trial);

        coefPP_current = zoeppritz(m_current(1), m_current(3), m_current(5), m_current(2), m_current(4), m_current(6), 1, 1, 0, anginc_PP);
        Rpp_current = real(coefPP_current);
        
        if forcepol
            cangle_trial = anginc_PP(find(coefPP_trial>0,1,'last'));
            cangle_current = anginc_PP(find(coefPP_current>0,1,'last'));
        end

    %calculate Rps for current model
        coefPS_trial = zoeppritz(m_trial(1), m_trial(3), m_trial(5), m_trial(2), m_trial(4), m_trial(6), 1, 2, 0, anginc_PS);
        Rps_trial = real(coefPS_trial); 

        coefPS_current = zoeppritz(m_current(1), m_current(3), m_current(5), m_current(2), m_current(4), m_current(6), 1, 2, 0, anginc_PS);
        Rps_current = real(coefPS_current);

    %calculate Rss for current model
        coefSS_trial = zoeppritz(m_trial(1), m_trial(3), m_trial(5), m_trial(2), m_trial(4), m_trial(6), 2, 2, 0, anginc_SS);
        Rss_trial = real(coefSS_trial); 

        coefSS_current = zoeppritz(m_current(1), m_current(3), m_current(5), m_current(2), m_current(4), m_current(6), 2, 2, 0, anginc_SS);
        Rss_current = real(coefSS_current);


    %%    likelihoods
    % calculate likelihoods

        if includePS == 1 && includeSS == 0

            %calculate weighted residuals based on a split normal
            %distribution
            
            r_t = [Rpp_input - Rpp_trial, Rps_input - Rps_trial]; %trial residuals
            rw_t = [(Rpp_input - Rpp_trial)./sigs_PP, (Rps_input - Rps_trial)./sigs_PS];

            r_c = [Rpp_input - Rpp_current, Rps_input - Rps_current];
            rw_c = [(Rpp_input - Rpp_current)./sigs_PP, (Rps_input - Rps_current)./sigs_PS];

            log_likelihood_trial =  - 0.5*nansum(rw_t.^2);
            log_likelihood_current = - 0.5*nansum(rw_c.^2);

        elseif includePS == 1 && includePP == 0
            r_t = [Rps_input - Rps_trial]; %trial residuals
            rw_t = [(Rps_input - Rps_trial)./sigs_PS];

            r_c = [Rps_input - Rps_current];
            rw_c = [(Rps_input - Rps_current)./sigs_PS];

            log_likelihood_trial =  - 0.5*nansum(rw_t.^2);
            log_likelihood_current = - 0.5*nansum(rw_c.^2);


        elseif includePS == 1 && includeSS == 1

            r_t = [Rpp_input - Rpp_trial, Rps_input - Rps_trial, Rss_input - Rss_trial];
            rw_t = [(Rpp_input - Rpp_trial)./sigs_PP, (Rps_input - Rps_trial)./sigs_PS, (Rss_input - Rss_trial)./sigs_SS];

            %log_likelihood_trial =  - 0.5*nansum(rw_t.^2);

            r_c = [Rpp_input - Rpp_current, Rps_input - Rps_current, Rss_input - Rss_current];
            rw_c = [(Rpp_input - Rpp_current)./sigs_PP, (Rps_input - Rps_current)./sigs_PS, (Rss_input - Rss_current)./sigs_SS];

            %log_likelihood_current = - 0.5*nansum(rw_c.^2);

            log_likelihood_trial =  - 0.5*nansum(rw_t.^2);
            log_likelihood_current = - 0.5*nansum(rw_c.^2);


        else
%{
            for jj=1:length(Rpp_input)
                if (Rpp_trial(jj)) > (Rpp_input(jj))
                    sigs_t_PP(jj) = sigsup_PP(jj);
                else
                    sigs_t_PP(jj) = sigsdn_PP(jj);
                end
            end
            for j=1:length(Rpp_input)
                if (Rpp_current(jj)) > (Rpp_input(jj))
                    sigs_c_PP(jj) = sigsup_PP(jj);
                else
                    sigs_c_PP(jj) = sigsdn_PP(jj);
                end
            end
%}
            r_t = Rpp_input - Rpp_trial; % trial residuals 
            rw_t = (Rpp_input - Rpp_trial)./sigs_PP;

            log_likelihood_trial =  - 0.5*nansum(rw_t.^2); % n.b. this does not include an error term yet

            r_c = Rpp_input - Rpp_current; % current residuals
            rw_c = (Rpp_input - Rpp_current)./sigs_PP;

            log_likelihood_current = - 0.5*nansum(rw_c.^2); % n.b. no error term yet

            log_likelihood_trial =  - 0.5*nansum(rw_t.^2);
            log_likelihood_current = - 0.5*nansum(rw_c.^2);
      

        end

    %% priors

    % calculate logs of the priors
        log_trial_prior_rho1 = log(1/(sqrt(2*pi)*std_rho1)) - ((m_trial(1)-rho1)^2)/(2*std_rho1^2);
        log_trial_prior_alpha1 = log(1/(sqrt(2*pi)*std_alpha1)) - ((m_trial(3)-alpha1).^2)/(2.*std_alpha1.^2);
        log_trial_prior_beta1 = log(1/(sqrt(2*pi)*std_beta1)) - ((m_trial(5)-beta1).^2)/(2.*std_beta1.^2);

        log_current_prior_rho1 = log(1/(sqrt(2*pi)*std_rho1)) - ((m_current(1)-rho1)^2)/(2*std_rho1^2);
        log_current_prior_alpha1 = log(1/(sqrt(2*pi)*std_alpha1)) - ((m_current(3)-alpha1).^2)/(2.*std_alpha1.^2);
        log_current_prior_beta1 = log(1/(sqrt(2*pi)*std_beta1)) - ((m_current(5)-beta1).^2)/(2.*std_beta1.^2);

    % combine priors of ice properties
        log_prior_trial = log_trial_prior_rho1 + log_trial_prior_alpha1 + log_trial_prior_beta1;
        log_prior_current = log_current_prior_rho1 + log_current_prior_alpha1 + log_current_prior_beta1;

        if forcepol
            log_trial_prior_cangle = log(1/(sqrt(2*pi)*cangle_std)) - ((cangle_trial-cangle_forced)^2)/(2*cangle_std^2);
            log_current_prior_cangle = log(1/(sqrt(2*pi)*cangle_std)) - ((cangle_current-cangle_forced)^2)/(2*cangle_std^2);

            log_prior_trial = log_prior_trial + log_trial_prior_cangle;
            log_prior_current = log_prior_current + log_current_prior_cangle;
        end

    %% posteriors
    %calc log posteriors

    log_posterior_trial = log_likelihood_trial + log_prior_trial;
    log_posterior_current = log_likelihood_current + log_prior_current;

    %% ratio of likelihoods
    % take ratio of likelihoods
        %log_ratio_likelihoods = log_likelihood_trial - log_current_priors_comb;


        log_ratio_posterior = log_posterior_trial - log_posterior_current;
        log_ratio_posterior_saved(i)=log_ratio_posterior;

        alphabeta_trial = m_trial(4)/m_trial(6);
        poisson_trial = (alphabeta_trial^2 - 2)/(2*(alphabeta_trial^2 - 1));
        % if model lies outside uniform priors for bed properties then reject
        if (m_trial(2) <= lims_rho2(1) || m_trial(2) >= lims_rho2(2) || m_trial(4) <= lims_alpha2(1) || m_trial(4) >= lims_alpha2(2) || m_trial(6) <= lims_beta2(1) || m_trial(6) >= lims_beta2(2))
            n_reject = n_reject + 1;
            accept_model = FALSE;
        elseif poisson_prior && (poisson_trial <= lims_poisson(1) || poisson_trial >= lims_poisson(2))
            n_reject = n_reject + 1;
            accept_model = FALSE;
        else
            if exp(log_ratio_posterior) > 1 %if t
                m_current = m_trial;
                n_accept = n_accept + 1;
                accept_model = TRUE;
            elseif exp(log_ratio_posterior)>rand
                n_accept = n_accept + 1;
                accept_model = TRUE;
            else
                n_reject = n_reject + 1;
                accept_model = FALSE;
            end

        end

        if accept_model
            m_current = m_trial;
            models_saved(i, :) = m_trial;
            log_likelihood_saved(i) = log_likelihood_trial;
            log_posterior_saved(i) = log_posterior_trial;
            log_prior_saved(i) = log_prior_trial;
            misfit_saved(i) = nansum(rw_t);
        else    
            models_saved(i, :) = m_current;
            log_likelihood_saved(i) = log_likelihood_current;
            log_posterior_saved(i) = log_posterior_current;
            log_prior_saved(i) = log_prior_current;
            misfit_saved(i) = nansum(rw_c);
        end

    end

    x = models_saved(burn_in+1:number_iterations, 4);
    y = models_saved(burn_in+1:number_iterations, 6);
    z = models_saved(burn_in+1:number_iterations, 2);
    %c = log_likelihood_saved(burn_in+1:number_iterations);
    c = exp(log_posterior_saved(burn_in+1:number_iterations));


    %[log_likelihood_best, besti] = max(log_likelihood_saved);
    [log_posterior_best, besti] = max(log_posterior_saved);
    m_best = models_saved(besti, :);

    %models_pbi = post burn in
    models_pbi = models_saved(burn_in+1:thinning:number_iterations, :);
    log_posterior_pbi = log_posterior_saved(burn_in+1:thinning:number_iterations);

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

    Z_saved = models_pbi(:,4).*models_pbi(:,2);
    alphabeta_saved = models_pbi(:,4)./models_pbi(:,6);
    poisson_saved = (alphabeta_saved.^2 - 2)./(2.*(alphabeta_saved.^2 - 1));

    %% find median model

    median_model = median(models_pbi);
    iqr_model = iqr(models_pbi);

    coefPP_med = zoeppritz(median_model(1), median_model(3), median_model(5), median_model(2), median_model(4), median_model(6), 1, 1, 0, anginc_PP);
    Rpp_med = real(coefPP_med);

    coefPS_med = zoeppritz(median_model(1), median_model(3), median_model(5), median_model(2), median_model(4), median_model(6), 1, 2, 0, anginc_PS);
    Rps_med = real(coefPS_med);

    coefSS_med = zoeppritz(median_model(1), median_model(3), median_model(5), median_model(2), median_model(4), median_model(6), 2, 2, 0, anginc_SS);
    Rss_med = real(coefSS_med);



    %% print out answers

    runningmedian = NaN(number_iterations-burn_in,6);
    for i=1:(number_iterations-burn_in)/1000
        runningmedianZ(i) = median(Z_saved(1:i*100));
        runningmedianpois(i) = median(poisson_saved(1:i*100));
        runningmedianbeta(i) = median(models_pbi(1:i*100, 6));
    end

    fprintf('\n----------------------------------------------------------------');
    if includeSS
        fprintf('\ntriple joint inversion');
    elseif includePS
        fprintf('\ndouble joint inversion');
    else
        fprintf('\nPP inversion');
    end
    fprintf('\n----------------------------------------------------------------');

    fprintf('\nbest model:');
    fprintf('\nVp2 = %.0f \x00B1 %.0f, Vs2 = %.0f \x00B1 %.0f, rho2 = %.0f \x00B1 %.0f', m_best(4), pm_of_models(4), m_best(6), pm_of_models(6), m_best(2), pm_of_models(2));
    fprintf('\nVp1 = %.0f \x00B1 %.0f, Vs1 = %.0f \x00B1 %.0f, rho1 = %.0f \x00B1 %.0f', m_best(3), pm_of_models(3), m_best(5), pm_of_models(5), m_best(1), pm_of_models(1));
    fprintf('\ncredible interval %.0f%%', cred_interval*100);
    fprintf('\nlog posterior %f', log_posterior_best);
    fprintf('\n----------------------------------------------------------------');

    fprintf('\n----------------------------------------------------------------');

    fprintf('\nmedian \x00B1 quartiles:');
    fprintf('\nVp2 = %.0f \x00B1 %.0f, Vs2 = %.0f \x00B1 %.0f, rho2 = %.0f \x00B1 %.0f', ...
        median(models_pbi(:,4)), 0.5*iqr(models_pbi(:,4)), median(models_pbi(:,6)), 0.5*iqr(models_pbi(:,6)), median(models_pbi(:,2)), 0.5*iqr(models_pbi(:,2)));
    fprintf('\nVp1 = %.0f \x00B1 %.0f, Vs1 = %.0f \x00B1 %.0f, rho1 = %.0f \x00B1 %.0f', median(models_pbi(:,3)), 0.5*iqr(models_pbi(:,3)), median(models_pbi(:,5)), 0.5*iqr(models_pbi(:,5)), median(models_pbi(:,1)), 0.5*iqr(models_pbi(:,1)));

    fprintf('\n----------------------------------------------------------------');
    
    fprintf('\n----------------------------------------------------------------');

    fprintf('\nmean \x00B1 std model:');
    fprintf('\nVp2 = %.0f \x00B1 %.0f, Vs2 = %.0f \x00B1 %.0f, rho2 = %.0f \x00B1 %.0f', ...
        mean(models_pbi(:,4)), std(models_pbi(:,4)), mean(models_pbi(:,6)), std(models_pbi(:,6)), mean(models_pbi(:,2)), std(models_pbi(:,2)));
    fprintf('\nVp1 = %.0f \x00B1 %.0f, Vs1 = %.0f \x00B1 %.0f, rho1 = %.0f \x00B1 %.0f', mean(models_pbi(:,3)), std(models_pbi(:,3)), mean(models_pbi(:,5)), std(models_pbi(:,5)), mean(models_pbi(:,1)), std(models_pbi(:,1)));

    fprintf('\n----------------------------------------------------------------');
    
    fprintf('\n----------------------------------------------------------------');

    
    fprintf('\nZ = (%.3f \x00B1 %.3f)e6, poisson = %.4f \x00B1 %.4f', median(Z_saved)*1e-6, 0.5*iqr(Z_saved)*1e-6, median(poisson_saved), 0.5*iqr(poisson_saved));
    fprintf('\nbest model: Z = %.3fe6, poisson = %.4f', Z_best*1e-6, poisson_best)
    fprintf('\nbest model: Z = %.3fe6, Vs = %.0f', Z_best*1e-6, m_best(6))
    fprintf('\nmean model: Z = %.3fe6, Vs = %.0f', Z_best*1e-6, mean(models_pbi(:,6)))
    fprintf('\n');
    fprintf('\nmodels accepted: %d, models rejected: %d', n_accept, n_reject);
    fprintf('\nfraction accepted %.2f',n_accept/number_iterations);
    fprintf('\n');

    %% plot results: posteriors, running medians, in terms of Z/poisson, and vp/vs/rho, and AVA of best fitting models
    plot_all_results_v3

    %% save results

    if (save_data == 1)&& (n_ensemble == 1)

    save_inversion_results;

    m_best_ensemble(k,:) = m_best;

    end
end
%% save if more than one sim in ensemble
if n_ensemble > 1
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
       
    priors.std_rho1 = std_rho1;
    priors.std_alpha1 = std_alpha1;
    priors.std_beta1 = std_beta1;
    priors.lims_rho2 = lims_rho2;
    priors.lims_alpha2 = lims_alpha2;
    priors.lims_beta2 = lims_beta2;
    priors.lims_poisson = lims_poisson;
    
    time = datestr(now, 'yyyy-mm-dd-hh-MM');
    formatOut = 'yyyy-mm-dd';
    foldername = strcat('inversion_results/ensemble/',datestr(now, formatOut));
    savedatafile = strcat(foldername,'/inversion_',time,'.mat');
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    
    
    
    ensemble_mean = mean(m_best_ensemble);
    ensemble_std = std(m_best_ensemble);
    alphabeta_ensemble = m_best_ensemble(:,4)./m_best_ensemble(:,6);
    poisson_ensemble = (alphabeta_ensemble.^2 - 2)./(2.*(alphabeta_ensemble.^2 - 1));
    Z_ensemble = m_best_ensemble(:,4).*m_best_ensemble(:,2);
    Z_ensemble_mean = mean(Z_ensemble);
    poisson_ensemble_mean = mean(poisson_ensemble);
    Z_ensemble_std = std(Z_ensemble);
    poisson_ensemble_std = std(poisson_ensemble);
    
    %%save it
    save(savedatafile, 'inputs', 'm_best_ensemble', 'ensemble_mean', 'ensemble_std', 'poisson_ensemble', 'Z_ensemble', 'priors');
    fprintf('\nsaving into file\n%s\n\n', savedatafile);
    
    
    %%plot Z and poisson for the ensemble mean
    figure; hold on;
    plot_zpois_misfit(Z_saved, poisson_saved, log_posterior_pbi)

    plot_zpois_ref(Z_ensemble_mean, poisson_ensemble_mean, poisson_ensemble_std, Z_ensemble_std)
    %plot

    coefPP_ens = zoeppritz(ensemble_mean(1), ensemble_mean(3), ensemble_mean(5), ensemble_mean(2), ensemble_mean(4), ensemble_mean(6), 1, 1, 0, anginc_PP);
    Rpp_ens = real(coefPP_ens);


    coefPS_ens = zoeppritz(ensemble_mean(1), ensemble_mean(3), ensemble_mean(5), ensemble_mean(2), ensemble_mean(4), ensemble_mean(6), 1, 2, 0, anginc_PS);
    Rps_ens = real(coefPS_ens); 
    
    %plot avo curves
    figure;
    subplot(1,2,1);
    hold on;
    plot(anginc_PP, Rpp_input, '.', 'Color','b');
    plot(anginc_PP, Rpp_ens, 'Color','b', 'LineWidth', 1.5);
    legend('data', 'ensemble solution');
    xlabel('ang inc (deg)'); ylabel('Rpp')

    subplot(1,2,2);
    hold on;
    plot(anginc_PS, Rps_input,'.', 'Color','r');
    plot(anginc_PS, Rps_ens, 'Color','r', 'LineWidth', 1.5);
    legend('data', 'ensemble solution');
    xlabel('ang inc (deg)'); ylabel('Rps')

    set(gcf, 'Name', '1_AVO_curves', 'Color','w');

end 
toc