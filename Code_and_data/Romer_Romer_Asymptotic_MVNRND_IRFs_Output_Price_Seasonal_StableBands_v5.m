%% Four Cholesky VAR specifications with seasonal dummies and stable-draw confidence bands
%
% Model suffixes:
%   _newm      = Romer-Romer new measure, no commodity prices
%   _newm_com  = Romer-Romer new measure, with commodity prices
%   _ff        = Federal funds rate, no commodity prices
%   _ff_com    = Federal funds rate, with commodity prices
%
% The OLS Cholesky IRF is the point estimate. Dashed lines show:
%
%   OLS IRF +/- one pointwise simulated standard error.
%
% The script also stores 95% pointwise normal-approximation confidence
% intervals for the requested responses and horizons. All models return and
% plot only the output and price-level responses. Policy variables and
% commodity prices remain in the estimated VARs but their responses are not
% returned, stored or plotted. Eleven monthly seasonal dummies are included
% as deterministic regressors in every equation of every VAR.
% The original OLS IRF is always retained as the point estimate, even if
% the OLS VAR is dynamically unstable. Confidence bands are computed only
% from stable asymptotic parameter draws with companion roots below one.

clear;
clc;
close all;

rng(4201, 'twister');

%% Requirements
%
% Place var_irf_asymptotic_mvnrnd_seasonal_stable_bands_v5.m in the same folder as this script.
% Statistics and Machine Learning Toolbox is required for MVNRND/WISHRND.
% The lag-selection section additionally requires Econometrics Toolbox.

%% Import the Romer and Romer data

data_table = readtable( ...
    "RomerandRomerDataAppendix.xls", ...
    'Sheet', "DATA BY MONTH");

disp(data_table.Properties.VariableNames);

%% Extract common series

industrial_production_log = data_table.LNIPNSA;
producer_price_log = data_table.LNPPINSA;
commodity_prices_log = data_table.LNWCP;
new_measure = data_table.SUMSHCK;
fed_funds = data_table.FF;

same_shock_series = isequal(data_table.SUMSHCK, data_table.SUMSHCKF);
fprintf('Are SUMSHCK and SUMSHCKF identical? %d\n', same_shock_series);

%% Construct monthly seasonal dummies
%
% An intercept and all 12 monthly dummies would be perfectly collinear.
% Therefore, the first seasonal category is omitted and 11 dummies are
% included in every VAR equation. When a valid month-number or datetime
% variable is available in the imported table, calendar months are used.
% Otherwise, the code uses the sequential monthly row order.

[seasonal_dummies, seasonal_dummy_names, month_of_year, ...
    seasonal_dummy_source] = build_monthly_seasonal_dummies(data_table);

fprintf('Seasonal-dummy source: %s\n', seasonal_dummy_source);
fprintf('Number of seasonal dummies included in each equation: %d\n', ...
    size(seasonal_dummies, 2));

%% Common inference settings

number_of_parameter_draws = 500;
maximum_irf_horizon = 48;
desired_policy_impact = 1;
confidence_level = 0.95;

% true draws both B and Sigma. false draws B with MVNRND and holds Sigma fixed.
draw_covariance_matrix = true;

% Stability condition for the draws used to construct the confidence bands.
% A draw is retained only when its largest absolute companion root is
% strictly below this threshold. Setting the threshold to 1.0 excludes all
% dynamically unstable draws while keeping the original OLS IRF unchanged.
stability_threshold = 1.0;

% Secondary numerical safeguard. This is applied only after the stability
% check and removes a retained stable draw if its normalised finite-horizon
% IRF is numerically extreme in the original model units.
maximum_allowed_absolute_irf = 100;

% Continue drawing until the requested number of stable, valid draws is
% retained. Increase this limit if the stable-draw acceptance rate is low.
maximum_draw_attempts = 500000;
draw_batch_size = 500;

asymptotic_options = struct( ...
    'desired_policy_impact', desired_policy_impact, ...
    'draw_covariance_matrix', draw_covariance_matrix, ...
    'confidence_level', confidence_level, ...
    'deterministic_regressors', seasonal_dummies, ...
    'deterministic_names', seasonal_dummy_names, ...
    'stability_threshold', stability_threshold, ...
    'maximum_allowed_absolute_irf', maximum_allowed_absolute_irf, ...
    'maximum_draw_attempts', maximum_draw_attempts, ...
    'draw_batch_size', draw_batch_size);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. New monetary policy measure VAR without commodity prices: _newm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The VAR includes an intercept and 11 monthly seasonal dummies.
% Cholesky ordering:
% 1. Log industrial production
% 2. Log producer price level
% 3. Cumulated Romer-Romer new monetary policy measure

y_newm = [ ...
    industrial_production_log, ...
    producer_price_log, ...
    new_measure];

lags_newm = 36;
policy_variable_index_newm = 3;

% Return only output and the price level. The policy response is still used
% internally to normalise the structural policy shock on impact.
asymptotic_options_newm = asymptotic_options;
asymptotic_options_newm.response_indices = [1, 2];

VAR_newm = var_irf_asymptotic_mvnrnd_seasonal_stable_bands_v5( ...
    y_newm, lags_newm, maximum_irf_horizon, ...
    policy_variable_index_newm, number_of_parameter_draws, ...
    asymptotic_options_newm);

horizons_newm = VAR_newm.horizons;
irfs_newm = VAR_newm.irf_ols;
irf_draws_newm = VAR_newm.irf_draws;
irfs_mean_draws_newm = VAR_newm.irf_mean_across_draws;
irfs_pointwise_std_newm = VAR_newm.irf_pointwise_std;
irfs_lower_1se_newm = VAR_newm.irf_lower_1se;
irfs_upper_1se_newm = VAR_newm.irf_upper_1se;
irfs_lower_95_newm = VAR_newm.irf_lower_normal;
irfs_upper_95_newm = VAR_newm.irf_upper_normal;

plot_scale_newm = [100; 100];

irfs_for_plot_newm = irfs_newm .* plot_scale_newm;
irfs_pointwise_std_plot_newm = ...
    irfs_pointwise_std_newm .* plot_scale_newm;
irfs_lower_1se_plot_newm = irfs_lower_1se_newm .* plot_scale_newm;
irfs_upper_1se_plot_newm = irfs_upper_1se_newm .* plot_scale_newm;
irfs_lower_95_plot_newm = irfs_lower_95_newm .* plot_scale_newm;
irfs_upper_95_plot_newm = irfs_upper_95_newm .* plot_scale_newm;

pointwise_std_table_newm = array2table( ...
    [horizons_newm', irfs_pointwise_std_plot_newm'], ...
    'VariableNames', {'Horizon', 'Output', 'PriceLevel'});

fprintf('\n_newm original policy impact response: %.6f\n', ...
    VAR_newm.policy_impact_before_normalisation);
fprintf('_newm OLS normalisation factor: %.6f\n', ...
    VAR_newm.normalisation_scaling_factor);
fprintf('_newm retained parameter draws: %d of %d\n', ...
    VAR_newm.accepted_number_of_draws, ...
    VAR_newm.requested_number_of_draws);
fprintf('_newm OLS maximum companion root: %.6f\n', ...
    VAR_newm.maximum_companion_root_ols);
fprintf('_newm stable-draw threshold for bands: %.6f\n', ...
    VAR_newm.stability_threshold);
fprintf(['_newm draw attempts: %d; rejected unstable draws: %d; ', ...
    'rejected extreme IRFs: %d; rejected invalid: %d; ', ...
    'acceptance rate: %.2f%%\n'], ...
    VAR_newm.attempted_number_of_draws, ...
    VAR_newm.rejected_unstable_draws, ...
    VAR_newm.rejected_extreme_irf_draws, ...
    VAR_newm.rejected_invalid_draws, ...
    100 * VAR_newm.draw_acceptance_rate);

% Returned response rows: 1 = output, 2 = price level.
panel_indices_newm = [1, 2];
panel_titles_newm = { ...
    'Effect on Output'; ...
    'Effect on the Price Level'};
y_axis_labels_newm = {'Percent'; 'Percent'};
y_axis_limits_newm = {[-5, 2]; [-7, 1]};
y_axis_ticks_newm = {-5:1:2; -7:1:1};

% Whole-figure title: edit this text to change the chart title.
figure_title_newm = ...
    'Romer-Romer New Monetary Policy Measure VAR';

[figure_newm, plot_layout_newm, axes_newm] = ...
    plot_tagged_irfs_with_bands( ...
    irfs_for_plot_newm, ...
    irfs_lower_1se_plot_newm, ...
    irfs_upper_1se_plot_newm, ...
    horizons_newm, ...
    panel_indices_newm, ...
    panel_titles_newm, ...
    y_axis_labels_newm, ...
    y_axis_limits_newm, ...
    y_axis_ticks_newm, ...
    [100, 50, 750, 750], ...
    figure_title_newm);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. New monetary policy measure VAR with commodity prices: _newm_com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The VAR includes an intercept and 11 monthly seasonal dummies.
% Cholesky ordering:
% 1. Log industrial production
% 2. Log producer price level
% 3. Log world commodity prices
% 4. Cumulated Romer-Romer new monetary policy measure

y_newm_com = [ ...
    industrial_production_log, ...
    producer_price_log, ...
    commodity_prices_log, ...
    new_measure];

lags_newm_com = 36;
policy_variable_index_newm_com = 4;

% Commodity prices and the policy variable remain in the VAR, but only
% output and the price level are returned. The policy response is still used
% internally to normalise the structural policy shock on impact.
asymptotic_options_newm_com = asymptotic_options;
asymptotic_options_newm_com.response_indices = [1, 2];

VAR_newm_com = var_irf_asymptotic_mvnrnd_seasonal_stable_bands_v5( ...
    y_newm_com, lags_newm_com, maximum_irf_horizon, ...
    policy_variable_index_newm_com, number_of_parameter_draws, ...
    asymptotic_options_newm_com);

horizons_newm_com = VAR_newm_com.horizons;
% Rows below are output and the price level only.
irfs_newm_com = VAR_newm_com.irf_ols;
irf_draws_newm_com = VAR_newm_com.irf_draws;
irfs_mean_draws_newm_com = VAR_newm_com.irf_mean_across_draws;
irfs_pointwise_std_newm_com = VAR_newm_com.irf_pointwise_std;
irfs_lower_1se_newm_com = VAR_newm_com.irf_lower_1se;
irfs_upper_1se_newm_com = VAR_newm_com.irf_upper_1se;
irfs_lower_95_newm_com = VAR_newm_com.irf_lower_normal;
irfs_upper_95_newm_com = VAR_newm_com.irf_upper_normal;

plot_scale_newm_com = [100; 100];

irfs_for_plot_newm_com = irfs_newm_com .* plot_scale_newm_com;
irfs_pointwise_std_plot_newm_com = ...
    irfs_pointwise_std_newm_com .* plot_scale_newm_com;
irfs_lower_1se_plot_newm_com = ...
    irfs_lower_1se_newm_com .* plot_scale_newm_com;
irfs_upper_1se_plot_newm_com = ...
    irfs_upper_1se_newm_com .* plot_scale_newm_com;
irfs_lower_95_plot_newm_com = ...
    irfs_lower_95_newm_com .* plot_scale_newm_com;
irfs_upper_95_plot_newm_com = ...
    irfs_upper_95_newm_com .* plot_scale_newm_com;

pointwise_std_table_newm_com = array2table( ...
    [horizons_newm_com', irfs_pointwise_std_plot_newm_com'], ...
    'VariableNames', {'Horizon', 'Output', 'PriceLevel'});

fprintf('\n_newm_com original policy impact response: %.6f\n', ...
    VAR_newm_com.policy_impact_before_normalisation);
fprintf('_newm_com OLS normalisation factor: %.6f\n', ...
    VAR_newm_com.normalisation_scaling_factor);
fprintf('_newm_com retained parameter draws: %d of %d\n', ...
    VAR_newm_com.accepted_number_of_draws, ...
    VAR_newm_com.requested_number_of_draws);
fprintf('_newm_com OLS maximum companion root: %.6f\n', ...
    VAR_newm_com.maximum_companion_root_ols);
fprintf('_newm_com stable-draw threshold for bands: %.6f\n', ...
    VAR_newm_com.stability_threshold);
fprintf(['_newm_com draw attempts: %d; rejected unstable draws: %d; ', ...
    'rejected extreme IRFs: %d; rejected invalid: %d; ', ...
    'acceptance rate: %.2f%%\n'], ...
    VAR_newm_com.attempted_number_of_draws, ...
    VAR_newm_com.rejected_unstable_draws, ...
    VAR_newm_com.rejected_extreme_irf_draws, ...
    VAR_newm_com.rejected_invalid_draws, ...
    100 * VAR_newm_com.draw_acceptance_rate);

% Returned response rows: 1 = output, 2 = price level.
panel_indices_newm_com = [1, 2];
panel_titles_newm_com = { ...
    'Effect on Output'; ...
    'Effect on the Price Level'};
y_axis_labels_newm_com = {'Percent'; 'Percent'};
y_axis_limits_newm_com = {[-5, 2]; [-7, 1]};
y_axis_ticks_newm_com = {-5:1:2; -7:1:1};

% Whole-figure title: edit this text to change the chart title.
figure_title_newm_com = ...
    'Romer-Romer New Monetary Policy Measure VAR with Commodity Prices';

[figure_newm_com, plot_layout_newm_com, axes_newm_com] = ...
    plot_tagged_irfs_with_bands( ...
    irfs_for_plot_newm_com, ...
    irfs_lower_1se_plot_newm_com, ...
    irfs_upper_1se_plot_newm_com, ...
    horizons_newm_com, ...
    panel_indices_newm_com, ...
    panel_titles_newm_com, ...
    y_axis_labels_newm_com, ...
    y_axis_limits_newm_com, ...
    y_axis_ticks_newm_com, ...
    [900, 50, 800, 750], ...
    figure_title_newm_com);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Federal funds rate VAR without commodity prices: _ff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The VAR includes an intercept and 11 monthly seasonal dummies.
% Cholesky ordering:
% 1. Log industrial production
% 2. Log producer price level
% 3. Federal funds rate

y_ff = [ ...
    industrial_production_log, ...
    producer_price_log, ...
    fed_funds];

lags_ff = 36;
policy_variable_index_ff = 3;

% Return only output and the price level. The funds-rate response is still
% used internally to normalise the structural policy shock on impact.
asymptotic_options_ff = asymptotic_options;
asymptotic_options_ff.response_indices = [1, 2];

VAR_ff = var_irf_asymptotic_mvnrnd_seasonal_stable_bands_v5( ...
    y_ff, lags_ff, maximum_irf_horizon, ...
    policy_variable_index_ff, number_of_parameter_draws, ...
    asymptotic_options_ff);

horizons_ff = VAR_ff.horizons;
irfs_ff = VAR_ff.irf_ols;
irf_draws_ff = VAR_ff.irf_draws;
irfs_mean_draws_ff = VAR_ff.irf_mean_across_draws;
irfs_pointwise_std_ff = VAR_ff.irf_pointwise_std;
irfs_lower_1se_ff = VAR_ff.irf_lower_1se;
irfs_upper_1se_ff = VAR_ff.irf_upper_1se;
irfs_lower_95_ff = VAR_ff.irf_lower_normal;
irfs_upper_95_ff = VAR_ff.irf_upper_normal;

plot_scale_ff = [100; 100];

irfs_for_plot_ff = irfs_ff .* plot_scale_ff;
irfs_pointwise_std_plot_ff = irfs_pointwise_std_ff .* plot_scale_ff;
irfs_lower_1se_plot_ff = irfs_lower_1se_ff .* plot_scale_ff;
irfs_upper_1se_plot_ff = irfs_upper_1se_ff .* plot_scale_ff;
irfs_lower_95_plot_ff = irfs_lower_95_ff .* plot_scale_ff;
irfs_upper_95_plot_ff = irfs_upper_95_ff .* plot_scale_ff;

pointwise_std_table_ff = array2table( ...
    [horizons_ff', irfs_pointwise_std_plot_ff'], ...
    'VariableNames', {'Horizon', 'Output', 'PriceLevel'});

fprintf('\n_ff original policy impact response: %.6f\n', ...
    VAR_ff.policy_impact_before_normalisation);
fprintf('_ff OLS normalisation factor: %.6f\n', ...
    VAR_ff.normalisation_scaling_factor);
fprintf('_ff retained parameter draws: %d of %d\n', ...
    VAR_ff.accepted_number_of_draws, VAR_ff.requested_number_of_draws);
fprintf('_ff OLS maximum companion root: %.6f\n', ...
    VAR_ff.maximum_companion_root_ols);
fprintf('_ff stable-draw threshold for bands: %.6f\n', ...
    VAR_ff.stability_threshold);
fprintf(['_ff draw attempts: %d; rejected unstable draws: %d; ', ...
    'rejected extreme IRFs: %d; rejected invalid: %d; ', ...
    'acceptance rate: %.2f%%\n'], ...
    VAR_ff.attempted_number_of_draws, ...
    VAR_ff.rejected_unstable_draws, ...
    VAR_ff.rejected_extreme_irf_draws, ...
    VAR_ff.rejected_invalid_draws, ...
    100 * VAR_ff.draw_acceptance_rate);

% Returned response rows: 1 = output, 2 = price level.
panel_indices_ff = [1, 2];
panel_titles_ff = { ...
    'Effect on Output'; ...
    'Effect on the Price Level'};
y_axis_labels_ff = {'Percent'; 'Percent'};
y_axis_limits_ff = {[-5, 2]; [-7, 1]};
y_axis_ticks_ff = {-5:1:2; -7:1:1};

% Whole-figure title: edit this text to change the chart title.
figure_title_ff = ...
    'Federal Funds Rate VAR';

[figure_ff, plot_layout_ff, axes_ff] = ...
    plot_tagged_irfs_with_bands( ...
    irfs_for_plot_ff, ...
    irfs_lower_1se_plot_ff, ...
    irfs_upper_1se_plot_ff, ...
    horizons_ff, ...
    panel_indices_ff, ...
    panel_titles_ff, ...
    y_axis_labels_ff, ...
    y_axis_limits_ff, ...
    y_axis_ticks_ff, ...
    [100, 100, 750, 750], ...
    figure_title_ff);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Federal funds rate VAR with commodity prices: _ff_com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The VAR includes an intercept and 11 monthly seasonal dummies.
% Cholesky ordering:
% 1. Log industrial production
% 2. Log producer price level
% 3. Log world commodity prices
% 4. Federal funds rate

y_ff_com = [ ...
    industrial_production_log, ...
    producer_price_log, ...
    commodity_prices_log, ...
    fed_funds];

lags_ff_com = 36;
policy_variable_index_ff_com = 4;

% Commodity prices and the policy variable remain in the VAR, but only
% output and the price level are returned. The funds-rate response is still
% used internally to normalise the structural policy shock on impact.
asymptotic_options_ff_com = asymptotic_options;
asymptotic_options_ff_com.response_indices = [1, 2];

VAR_ff_com = var_irf_asymptotic_mvnrnd_seasonal_stable_bands_v5( ...
    y_ff_com, lags_ff_com, maximum_irf_horizon, ...
    policy_variable_index_ff_com, number_of_parameter_draws, ...
    asymptotic_options_ff_com);

horizons_ff_com = VAR_ff_com.horizons;
% Rows below are output and the price level only.
irfs_ff_com = VAR_ff_com.irf_ols;
irf_draws_ff_com = VAR_ff_com.irf_draws;
irfs_mean_draws_ff_com = VAR_ff_com.irf_mean_across_draws;
irfs_pointwise_std_ff_com = VAR_ff_com.irf_pointwise_std;
irfs_lower_1se_ff_com = VAR_ff_com.irf_lower_1se;
irfs_upper_1se_ff_com = VAR_ff_com.irf_upper_1se;
irfs_lower_95_ff_com = VAR_ff_com.irf_lower_normal;
irfs_upper_95_ff_com = VAR_ff_com.irf_upper_normal;

plot_scale_ff_com = [100; 100];

irfs_for_plot_ff_com = irfs_ff_com .* plot_scale_ff_com;
irfs_pointwise_std_plot_ff_com = ...
    irfs_pointwise_std_ff_com .* plot_scale_ff_com;
irfs_lower_1se_plot_ff_com = ...
    irfs_lower_1se_ff_com .* plot_scale_ff_com;
irfs_upper_1se_plot_ff_com = ...
    irfs_upper_1se_ff_com .* plot_scale_ff_com;
irfs_lower_95_plot_ff_com = ...
    irfs_lower_95_ff_com .* plot_scale_ff_com;
irfs_upper_95_plot_ff_com = ...
    irfs_upper_95_ff_com .* plot_scale_ff_com;

pointwise_std_table_ff_com = array2table( ...
    [horizons_ff_com', irfs_pointwise_std_plot_ff_com'], ...
    'VariableNames', {'Horizon', 'Output', 'PriceLevel'});

fprintf('\n_ff_com original policy impact response: %.6f\n', ...
    VAR_ff_com.policy_impact_before_normalisation);
fprintf('_ff_com OLS normalisation factor: %.6f\n', ...
    VAR_ff_com.normalisation_scaling_factor);
fprintf('_ff_com retained parameter draws: %d of %d\n', ...
    VAR_ff_com.accepted_number_of_draws, ...
    VAR_ff_com.requested_number_of_draws);
fprintf('_ff_com OLS maximum companion root: %.6f\n', ...
    VAR_ff_com.maximum_companion_root_ols);
fprintf('_ff_com stable-draw threshold for bands: %.6f\n', ...
    VAR_ff_com.stability_threshold);
fprintf(['_ff_com draw attempts: %d; rejected unstable draws: %d; ', ...
    'rejected extreme IRFs: %d; rejected invalid: %d; ', ...
    'acceptance rate: %.2f%%\n'], ...
    VAR_ff_com.attempted_number_of_draws, ...
    VAR_ff_com.rejected_unstable_draws, ...
    VAR_ff_com.rejected_extreme_irf_draws, ...
    VAR_ff_com.rejected_invalid_draws, ...
    100 * VAR_ff_com.draw_acceptance_rate);

% Returned response rows: 1 = output, 2 = price level.
panel_indices_ff_com = [1, 2];
panel_titles_ff_com = { ...
    'Effect on Output'; ...
    'Effect on the Price Level'};
y_axis_labels_ff_com = {'Percent'; 'Percent'};
y_axis_limits_ff_com = {[-5, 2]; [-7, 1]};
y_axis_ticks_ff_com = {-5:1:2; -7:1:1};

% Whole-figure title: edit this text to change the chart title.
figure_title_ff_com = ...
    'Federal Funds Rate VAR with Commodity Prices';

[figure_ff_com, plot_layout_ff_com, axes_ff_com] = ...
    plot_tagged_irfs_with_bands( ...
    irfs_for_plot_ff_com, ...
    irfs_lower_1se_plot_ff_com, ...
    irfs_upper_1se_plot_ff_com, ...
    horizons_ff_com, ...
    panel_indices_ff_com, ...
    panel_titles_ff_com, ...
    y_axis_labels_ff_com, ...
    y_axis_limits_ff_com, ...
    y_axis_ticks_ff_com, ...
    [900, 100, 800, 750], ...
    figure_title_ff_com);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assumption tests: VAR lag-length selection with seasonal dummies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y = [ ...
    industrial_production_log, ...
    producer_price_log, ...
    new_measure];

maximum_lag = 60;
alpha = 0.05;

if any(~isfinite(Y), 'all')
    error(['Y contains NaN or Inf values. Handle missing observations ', ...
        'before selecting the VAR lag length.']);
end

[number_of_observations, number_of_variables] = size(Y);
Y_presample = Y(1:maximum_lag, :);
Y_estimation = Y(maximum_lag + 1:end, :);
seasonal_dummies_estimation = ...
    seasonal_dummies(maximum_lag + 1:end, :);
effective_sample_size = size(Y_estimation, 1);

number_of_seasonal_dummies = size(seasonal_dummies, 2);
maximum_number_of_regressors = 1 + number_of_seasonal_dummies + ...
    number_of_variables * maximum_lag;

if effective_sample_size <= maximum_number_of_regressors
    error(['The maximum lag is too large relative to the sample size. ', ...
        'Reduce maximum_lag.']);
end

lag_orders = (1:maximum_lag)';
log_likelihood = NaN(maximum_lag, 1);
number_of_parameters = NaN(maximum_lag, 1);
estimated_models = cell(maximum_lag, 1);

for lag = 1:maximum_lag
    model = varm(number_of_variables, lag);
    [estimated_models{lag}, ~, log_likelihood(lag)] = estimate( ...
        model, Y_estimation, ...
        'Y0', Y_presample, ...
        'X', seasonal_dummies_estimation, ...
        'Display', 'off');
    model_summary = summarize(estimated_models{lag});
    number_of_parameters(lag) = model_summary.NumEstimatedParameters;
end

[AIC, SIC] = aicbic( ...
    log_likelihood, number_of_parameters, effective_sample_size);

HQC = -2 .* log_likelihood + ...
    2 .* number_of_parameters .* log(log(effective_sample_size));

AIC = AIC(:);
HQC = HQC(:);
SIC = SIC(:);

information_criteria_table = table( ...
    lag_orders, log_likelihood, number_of_parameters, AIC, HQC, SIC, ...
    'VariableNames', {'Lag', 'LogLikelihood', 'NumberOfParameters', ...
    'AIC', 'HQC', 'SIC_BIC'});

disp('Information criteria results:')
disp(information_criteria_table)

[~, AIC_selected_lag] = min(AIC);
[~, HQC_selected_lag] = min(HQC);
[~, SIC_selected_lag] = min(SIC);

fprintf('\nInformation-criteria lag selection:\n');
fprintf('AIC selects lag: %d\n', AIC_selected_lag);
fprintf('HQC selects lag: %d\n', HQC_selected_lag);
fprintf('SIC/BIC selects lag: %d\n', SIC_selected_lag);

restricted_lag = (1:maximum_lag - 1)';
unrestricted_lag = (2:maximum_lag)';
LR_statistic = NaN(maximum_lag - 1, 1);
LR_degrees_of_freedom = NaN(maximum_lag - 1, 1);
LR_p_value = NaN(maximum_lag - 1, 1);
LR_critical_value = NaN(maximum_lag - 1, 1);
LR_reject_restricted_model = false(maximum_lag - 1, 1);

for lag = 2:maximum_lag
    current_test = lag - 1;
    unrestricted_log_likelihood = log_likelihood(lag);
    restricted_log_likelihood = log_likelihood(lag - 1);
    degrees_of_freedom = ...
        number_of_parameters(lag) - number_of_parameters(lag - 1);

    [reject_null, p_value, statistic, critical_value] = lratiotest( ...
        unrestricted_log_likelihood, restricted_log_likelihood, ...
        degrees_of_freedom, alpha);

    LR_statistic(current_test) = statistic;
    LR_degrees_of_freedom(current_test) = degrees_of_freedom;
    LR_p_value(current_test) = p_value;
    LR_critical_value(current_test) = critical_value;
    LR_reject_restricted_model(current_test) = reject_null;
end

LR_test_table = table( ...
    restricted_lag, unrestricted_lag, LR_statistic, ...
    LR_degrees_of_freedom, LR_p_value, LR_critical_value, ...
    LR_reject_restricted_model, ...
    'VariableNames', {'RestrictedLag', 'UnrestrictedLag', 'LRStatistic', ...
    'DegreesOfFreedom', 'PValue', 'CriticalValue', ...
    'RejectRestrictedModel'});

disp(' ')
disp('Sequential likelihood-ratio tests:')
disp(LR_test_table)

LR_selected_lag = maximum_lag;

for lag = maximum_lag:-1:2
    current_test = lag - 1;
    if LR_reject_restricted_model(current_test)
        LR_selected_lag = lag;
        break
    else
        LR_selected_lag = lag - 1;
    end
end

fprintf('\nSequential LR testing selects lag: %d\n', LR_selected_lag);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local plotting function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [figure_handle, layout_handle, axes_handles] = ...
        plot_tagged_irfs_with_bands( ...
        irfs_for_plot, irfs_lower_band, irfs_upper_band, horizons, ...
        panel_indices, panel_titles, y_axis_labels, y_axis_limits, ...
        y_axis_ticks, figure_position, figure_title)
%PLOT_TAGGED_IRFS_WITH_BANDS Plot OLS IRFs and pointwise +/-1 SE bands.
% figure_title controls the editable title above the complete tiled chart.

    number_of_panels = numel(panel_indices);
    panel_letters = 'abcdefghijklmnopqrstuvwxyz';

    figure_handle = figure('Position', figure_position, 'Color', 'w');
    layout_handle = tiledlayout(number_of_panels, 1);
    layout_handle.TileSpacing = 'compact';
    layout_handle.Padding = 'compact';

    % Whole-chart title displayed above all impulse-response panels.
    title( ...
        layout_handle, ...
        figure_title, ...
        'FontSize', 20, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'none');

    axes_handles = gobjects(number_of_panels, 1);

    for panel = 1:number_of_panels
        axes_handles(panel) = nexttile(layout_handle);
        current_axis = axes_handles(panel);
        current_variable = panel_indices(panel);
        hold(current_axis, 'on');

        point_estimate_line = plot( ...
            current_axis, horizons, irfs_for_plot(current_variable, :), ...
            'k-', 'LineWidth', 1.5);

        lower_band_line = plot( ...
            current_axis, horizons, irfs_lower_band(current_variable, :), ...
            'k--', 'LineWidth', 1.1);

        plot( ...
            current_axis, horizons, irfs_upper_band(current_variable, :), ...
            'k--', 'LineWidth', 1.1);

        yline(current_axis, 0, 'k:', 'LineWidth', 1, ...
            'HandleVisibility', 'off');

        title(current_axis, sprintf('%s. %s', ...
            panel_letters(panel), panel_titles{panel}), 'FontSize', 18);
        xlabel(current_axis, 'Months after Shock', 'FontSize', 16);
        ylabel(current_axis, y_axis_labels{panel}, 'FontSize', 16);

        if ~isempty(y_axis_limits{panel})
            ylim(current_axis, y_axis_limits{panel});
        end

        if ~isempty(y_axis_ticks{panel})
            yticks(current_axis, y_axis_ticks{panel});
        end

        xlim(current_axis, [0, 48]);
        xticks(current_axis, 0:3:48);
        current_axis.FontSize = 14;
        current_axis.Box = 'on';
        current_axis.LineWidth = 1;

        if panel == 1
            legend(current_axis, [point_estimate_line, lower_band_line], ...
                {'OLS IRF', 'Pointwise \pm 1 SE'}, ...
                'Location', 'best', 'FontSize', 12);
        end

        hold(current_axis, 'off');
    end

    linkaxes(axes_handles, 'x');
end


function [seasonal_dummies, seasonal_dummy_names, month_of_year, source] = ...
        build_monthly_seasonal_dummies(data_table)
%BUILD_MONTHLY_SEASONAL_DUMMIES Create 11 monthly seasonal indicators.
%
% The function first searches for a numeric month variable with values 1-12.
% It then searches for a datetime variable. If neither is available, it uses
% the sequential row position, which is appropriate for a complete monthly
% data sheet. The first seasonal category is omitted to avoid collinearity
% with the intercept.

    number_of_observations = height(data_table);
    variable_names = data_table.Properties.VariableNames;
    lower_variable_names = lower(string(variable_names));

    month_of_year = [];
    source = '';

    candidate_month_names = [ ...
        "month", "mon", "monthnumber", "month_number", "mth"];

    for candidate = candidate_month_names
        location = find(lower_variable_names == candidate, 1);

        if isempty(location)
            continue
        end

        candidate_values = data_table.(variable_names{location});

        if isnumeric(candidate_values) && ...
                isvector(candidate_values) && ...
                numel(candidate_values) == number_of_observations && ...
                all(isfinite(candidate_values)) && ...
                all(candidate_values >= 1 & candidate_values <= 12) && ...
                all(candidate_values == round(candidate_values))

            month_of_year = double(candidate_values(:));
            source = sprintf('numeric table variable "%s"', ...
                variable_names{location});
            break
        end
    end

    if isempty(month_of_year)
        for variable = 1:numel(variable_names)
            candidate_values = data_table.(variable_names{variable});

            if isdatetime(candidate_values) && ...
                    isvector(candidate_values) && ...
                    numel(candidate_values) == number_of_observations && ...
                    all(~isnat(candidate_values))

                month_of_year = month(candidate_values(:));
                source = sprintf('datetime table variable "%s"', ...
                    variable_names{variable});
                break
            end
        end
    end

    if isempty(month_of_year)
        month_of_year = mod((0:number_of_observations - 1)', 12) + 1;
        source = ['sequential monthly row order; the omitted category is ' ...
            'the first month represented in the data'];
    end

    omitted_month = 1;
    included_months = setdiff(1:12, omitted_month, 'stable');
    seasonal_dummies = zeros(number_of_observations, 11);

    for dummy = 1:numel(included_months)
        seasonal_dummies(:, dummy) = ...
            double(month_of_year == included_months(dummy));
    end

    if startsWith(source, 'sequential monthly row order')
        seasonal_dummy_names = ...
            "SeasonalCategory" + string(included_months) + "Dummy";
    else
        calendar_month_names = [ ...
            "January", "February", "March", "April", "May", "June", ...
            "July", "August", "September", "October", "November", ...
            "December"];

        seasonal_dummy_names = ...
            calendar_month_names(included_months) + "Dummy";
    end

end

