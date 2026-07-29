%% Four Cholesky VAR specifications: output and price-level IRFs only
%
% This standalone script estimates four monthly Cholesky VARs:
%   1. Romer-Romer new measure, without commodity prices
%   2. Romer-Romer new measure, with commodity prices
%   3. Federal funds rate, without commodity prices
%   4. Federal funds rate, with commodity prices
%
% Every VAR includes an intercept and 11 monthly seasonal dummies.
% The policy and commodity-price variables remain in the estimated systems,
% but only the impulse responses of output and the price level are returned
% and plotted. No confidence bands or policy-variable responses are plotted.

clear;
clc;
close all;

%% Import the Romer and Romer data

data_table = readtable( ...
    "RomerandRomerDataAppendix.xls", ...
    'Sheet', "DATA BY MONTH");

disp(data_table.Properties.VariableNames);

%% Extract the series

industrial_production_log = data_table.LNIPNSA;
producer_price_log = data_table.LNPPINSA;
commodity_prices_log = data_table.LNWCP;
new_measure = data_table.SUMSHCK;
fed_funds = data_table.FF;

%% Construct 11 monthly seasonal dummies
%
% The first monthly category is omitted because the VAR also includes an
% intercept. This avoids perfect multicollinearity.

[seasonal_dummies, seasonal_dummy_names, seasonal_dummy_source] = ...
    build_monthly_seasonal_dummies(data_table);

fprintf('Seasonal-dummy source: %s\n', seasonal_dummy_source);
fprintf('Seasonal dummies included in each equation: %d\n', ...
    size(seasonal_dummies, 2));

disp('Seasonal dummy names:');
disp(seasonal_dummy_names');

%% Common settings

number_of_lags = 36;
maximum_irf_horizon = 48;
desired_policy_impact = 1;

% Only these responses are returned from each VAR:
%   1 = log industrial production
%   2 = log producer price level
response_indices = [1, 2];

% Both variables are in natural logarithms, so multiply their responses by
% 100 to express the log-point responses approximately as percentages.
plot_scale = [100; 100];

panel_titles = { ...
    'Effect on Output'; ...
    'Effect on the Price Level'};

y_axis_labels = {'Percent'; 'Percent'};
y_axis_limits = {[-5, 2]; [-7, 1]};
y_axis_ticks = {-5:1:2; -7:1:1};

%% 1. Romer-Romer new measure VAR without commodity prices
%
% Cholesky ordering:
%   1. Log industrial production
%   2. Log producer price level
%   3. Cumulated Romer-Romer monetary policy measure

y_newm = [ ...
    industrial_production_log, ...
    producer_price_log, ...
    new_measure];

policy_variable_index_newm = 3;

VAR_newm = estimate_cholesky_var_irfs( ...
    y_newm, ...
    number_of_lags, ...
    seasonal_dummies, ...
    policy_variable_index_newm, ...
    maximum_irf_horizon, ...
    desired_policy_impact, ...
    response_indices);

irfs_for_plot_newm = VAR_newm.irfs .* plot_scale;

fprintf('\n_newm policy impact before normalisation: %.6f\n', ...
    VAR_newm.policy_impact_before_normalisation);
fprintf('_newm normalisation factor: %.6f\n', ...
    VAR_newm.normalisation_factor);

[figure_newm, layout_newm, axes_newm] = plot_output_price_irfs( ...
    VAR_newm.horizons, ...
    irfs_for_plot_newm, ...
    panel_titles, ...
    y_axis_labels, ...
    y_axis_limits, ...
    y_axis_ticks, ...
    [100, 50, 750, 750], ...
    'Romer-Romer Measure VAR');

%% 2. Romer-Romer new measure VAR with commodity prices
%
% Cholesky ordering:
%   1. Log industrial production
%   2. Log producer price level
%   3. Log world commodity prices
%   4. Cumulated Romer-Romer monetary policy measure

y_newm_com = [ ...
    industrial_production_log, ...
    producer_price_log, ...
    commodity_prices_log, ...
    new_measure];

policy_variable_index_newm_com = 4;

VAR_newm_com = estimate_cholesky_var_irfs( ...
    y_newm_com, ...
    number_of_lags, ...
    seasonal_dummies, ...
    policy_variable_index_newm_com, ...
    maximum_irf_horizon, ...
    desired_policy_impact, ...
    response_indices);

irfs_for_plot_newm_com = VAR_newm_com.irfs .* plot_scale;

fprintf('\n_newm_com policy impact before normalisation: %.6f\n', ...
    VAR_newm_com.policy_impact_before_normalisation);
fprintf('_newm_com normalisation factor: %.6f\n', ...
    VAR_newm_com.normalisation_factor);

[figure_newm_com, layout_newm_com, axes_newm_com] = ...
    plot_output_price_irfs( ...
    VAR_newm_com.horizons, ...
    irfs_for_plot_newm_com, ...
    panel_titles, ...
    y_axis_labels, ...
    y_axis_limits, ...
    y_axis_ticks, ...
    [900, 50, 800, 750], ...
    'Romer-Romer Measure VAR with Commodity Prices');

%% 3. Federal funds rate VAR without commodity prices
%
% Cholesky ordering:
%   1. Log industrial production
%   2. Log producer price level
%   3. Federal funds rate

y_ff = [ ...
    industrial_production_log, ...
    producer_price_log, ...
    fed_funds];

policy_variable_index_ff = 3;

VAR_ff = estimate_cholesky_var_irfs( ...
    y_ff, ...
    number_of_lags, ...
    seasonal_dummies, ...
    policy_variable_index_ff, ...
    maximum_irf_horizon, ...
    desired_policy_impact, ...
    response_indices);

irfs_for_plot_ff = VAR_ff.irfs .* plot_scale;

fprintf('\n_ff policy impact before normalisation: %.6f\n', ...
    VAR_ff.policy_impact_before_normalisation);
fprintf('_ff normalisation factor: %.6f\n', ...
    VAR_ff.normalisation_factor);

[figure_ff, layout_ff, axes_ff] = plot_output_price_irfs( ...
    VAR_ff.horizons, ...
    irfs_for_plot_ff, ...
    panel_titles, ...
    y_axis_labels, ...
    y_axis_limits, ...
    y_axis_ticks, ...
    [100, 100, 750, 750], ...
    'Federal Funds Rate VAR');

%% 4. Federal funds rate VAR with commodity prices
%
% Cholesky ordering:
%   1. Log industrial production
%   2. Log producer price level
%   3. Log world commodity prices
%   4. Federal funds rate

y_ff_com = [ ...
    industrial_production_log, ...
    producer_price_log, ...
    commodity_prices_log, ...
    fed_funds];

policy_variable_index_ff_com = 4;

VAR_ff_com = estimate_cholesky_var_irfs( ...
    y_ff_com, ...
    number_of_lags, ...
    seasonal_dummies, ...
    policy_variable_index_ff_com, ...
    maximum_irf_horizon, ...
    desired_policy_impact, ...
    response_indices);

irfs_for_plot_ff_com = VAR_ff_com.irfs .* plot_scale;

fprintf('\n_ff_com policy impact before normalisation: %.6f\n', ...
    VAR_ff_com.policy_impact_before_normalisation);
fprintf('_ff_com normalisation factor: %.6f\n', ...
    VAR_ff_com.normalisation_factor);

[figure_ff_com, layout_ff_com, axes_ff_com] = plot_output_price_irfs( ...
    VAR_ff_com.horizons, ...
    irfs_for_plot_ff_com, ...
    panel_titles, ...
    y_axis_labels, ...
    y_axis_limits, ...
    y_axis_ticks, ...
    [900, 100, 800, 750], ...
    'Federal Funds Rate VAR with Commodity Prices');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local function: estimate the point-estimate Cholesky VAR IRFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = estimate_cholesky_var_irfs( ...
        Y, number_of_lags, deterministic_regressors, ...
        policy_variable_index, maximum_irf_horizon, ...
        desired_policy_impact, response_indices)
%ESTIMATE_CHOLESKY_VAR_IRFS Estimate a reduced-form VAR by OLS and calculate
%recursive Cholesky impulse responses to the selected policy innovation.
%
% Only the rows listed in response_indices are returned in results.irfs.

    if ~isnumeric(Y) || ndims(Y) ~= 2
        error('Y must be a numeric observations-by-variables matrix.');
    end

    [number_of_observations, number_of_variables] = size(Y);

    if number_of_lags < 1 || number_of_lags ~= floor(number_of_lags)
        error('number_of_lags must be a positive integer.');
    end

    if number_of_observations <= number_of_lags
        error('The sample must contain more observations than VAR lags.');
    end

    if size(deterministic_regressors, 1) ~= number_of_observations
        error(['The deterministic regressors must have the same number ', ...
            'of rows as Y.']);
    end

    if any(~isfinite(Y), 'all') || ...
            any(~isfinite(deterministic_regressors), 'all')
        error('Y and the deterministic regressors cannot contain NaN or Inf.');
    end

    if policy_variable_index < 1 || ...
            policy_variable_index > number_of_variables
        error('policy_variable_index is outside the columns of Y.');
    end

    if any(response_indices < 1) || ...
            any(response_indices > number_of_variables)
        error('response_indices contains an invalid variable index.');
    end

    effective_sample_size = number_of_observations - number_of_lags;
    number_of_deterministic_terms = ...
        1 + size(deterministic_regressors, 2);
    number_of_regressors = number_of_deterministic_terms + ...
        number_of_variables * number_of_lags;

    if effective_sample_size <= number_of_regressors
        error(['There are not enough usable observations for this VAR. ', ...
            'Reduce the number of lags or deterministic regressors.']);
    end

    dependent_variables = Y(number_of_lags + 1:end, :);

    % Intercept and contemporaneous deterministic regressors.
    regressors = [ ...
        ones(effective_sample_size, 1), ...
        deterministic_regressors(number_of_lags + 1:end, :)];

    % Add lagged values in blocks:
    % [Y(t-1), Y(t-2), ..., Y(t-p)].
    for lag = 1:number_of_lags
        lagged_block = Y( ...
            number_of_lags + 1 - lag:number_of_observations - lag, :);
        regressors = [regressors, lagged_block]; %#ok<AGROW>
    end

    % Each column of coefficient_matrix corresponds to one VAR equation.
    coefficient_matrix = regressors \ dependent_variables;
    residuals = dependent_variables - regressors * coefficient_matrix;

    residual_degrees_of_freedom = ...
        effective_sample_size - rank(regressors);

    if residual_degrees_of_freedom <= 0
        error('The estimated VAR has no positive residual degrees of freedom.');
    end

    residual_covariance = ...
        (residuals' * residuals) / residual_degrees_of_freedom;
    residual_covariance = ...
        (residual_covariance + residual_covariance') / 2;

    [impact_matrix, chol_flag] = chol(residual_covariance, 'lower');

    if chol_flag ~= 0
        error(['The residual covariance matrix is not positive definite. ', ...
            'Check the data, lag length and deterministic regressors.']);
    end

    % Recover the K-by-K autoregressive coefficient matrix for every lag.
    autoregressive_matrices = zeros( ...
        number_of_variables, number_of_variables, number_of_lags);

    first_lag_row = number_of_deterministic_terms + 1;

    for lag = 1:number_of_lags
        current_rows = first_lag_row + ...
            (lag - 1) * number_of_variables: ...
            first_lag_row + lag * number_of_variables - 1;

        autoregressive_matrices(:, :, lag) = ...
            coefficient_matrix(current_rows, :)';
    end

    % The relevant structural shock is the policy-variable column of the
    % lower-triangular Cholesky impact matrix.
    policy_impact_before_normalisation = ...
        impact_matrix(policy_variable_index, policy_variable_index);

    if abs(policy_impact_before_normalisation) < eps
        error('The policy innovation has a numerically zero impact response.');
    end

    normalisation_factor = ...
        desired_policy_impact / policy_impact_before_normalisation;

    impact_response = ...
        impact_matrix(:, policy_variable_index) * normalisation_factor;

    full_irfs = zeros(number_of_variables, maximum_irf_horizon + 1);
    full_irfs(:, 1) = impact_response;

    % VAR recursion:
    % response(h) = A1 response(h-1) + ... + Ap response(h-p).
    for horizon = 1:maximum_irf_horizon
        for lag = 1:min(number_of_lags, horizon)
            full_irfs(:, horizon + 1) = ...
                full_irfs(:, horizon + 1) + ...
                autoregressive_matrices(:, :, lag) * ...
                full_irfs(:, horizon - lag + 1);
        end
    end

    results = struct();
    results.horizons = 0:maximum_irf_horizon;
    results.irfs = full_irfs(response_indices, :);
    results.response_indices = response_indices;
    results.policy_impact_before_normalisation = ...
        policy_impact_before_normalisation;
    results.normalisation_factor = normalisation_factor;
    results.coefficient_matrix = coefficient_matrix;
    results.residuals = residuals;
    results.residual_covariance = residual_covariance;
    results.impact_matrix = impact_matrix;
    results.autoregressive_matrices = autoregressive_matrices;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local function: plot output and price-level responses only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [figure_handle, layout_handle, axes_handles] = ...
        plot_output_price_irfs( ...
        horizons, irfs_for_plot, panel_titles, y_axis_labels, ...
        y_axis_limits, y_axis_ticks, figure_position, figure_title)
%PLOT_OUTPUT_PRICE_IRFS Plot exactly two panels: output and price level.

    if size(irfs_for_plot, 1) ~= 2
        error(['irfs_for_plot must contain exactly two rows: ', ...
            'output followed by the price level.']);
    end

    number_of_panels = 2;
    panel_letters = 'ab';

    figure_handle = figure( ...
        'Position', figure_position, ...
        'Color', 'w');

    layout_handle = tiledlayout(number_of_panels, 1);
    layout_handle.TileSpacing = 'compact';
    layout_handle.Padding = 'compact';

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

        plot( ...
            current_axis, ...
            horizons, ...
            irfs_for_plot(panel, :), ...
            'k-', ...
            'LineWidth', 1.5);

        hold(current_axis, 'on');

        yline( ...
            current_axis, ...
            0, ...
            'k:', ...
            'LineWidth', 1, ...
            'HandleVisibility', 'off');

        hold(current_axis, 'off');

        title( ...
            current_axis, ...
            sprintf('%s. %s', panel_letters(panel), panel_titles{panel}), ...
            'FontSize', 18);

        xlabel( ...
            current_axis, ...
            'Months after Shock', ...
            'FontSize', 16);

        ylabel( ...
            current_axis, ...
            y_axis_labels{panel}, ...
            'FontSize', 16);

        if ~isempty(y_axis_limits{panel})
            ylim(current_axis, y_axis_limits{panel});
        end

        if ~isempty(y_axis_ticks{panel})
            yticks(current_axis, y_axis_ticks{panel});
        end

        xlim(current_axis, [0, max(horizons)]);
        xticks(current_axis, 0:3:max(horizons));

        current_axis.FontSize = 14;
        current_axis.Box = 'on';
        current_axis.LineWidth = 1;
    end

    linkaxes(axes_handles, 'x');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local function: construct monthly seasonal dummies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [seasonal_dummies, seasonal_dummy_names, source] = ...
        build_monthly_seasonal_dummies(data_table)
%BUILD_MONTHLY_SEASONAL_DUMMIES Create 11 monthly indicators.
%
% The function first searches for a numeric month variable containing values
% from 1 to 12. It then searches for a datetime variable. If neither exists,
% sequential monthly row order is used. Month/category 1 is omitted.

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
            source = sprintf( ...
                'numeric table variable "%s"', ...
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
                source = sprintf( ...
                    'datetime table variable "%s"', ...
                    variable_names{variable});
                break
            end
        end
    end

    if isempty(month_of_year)
        month_of_year = mod((0:number_of_observations - 1)', 12) + 1;
        source = ['sequential monthly row order; category 1 is omitted'];
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
