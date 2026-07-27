function results = var_irf_asymptotic_mvnrnd(Y, lags, maximum_horizon, policy_index, number_of_draws, options)
%VAR_IRF_ASYMPTOTIC_MVNRND Frequentist Monte Carlo bands for a Cholesky VAR.
%
% Estimates a reduced-form VAR by OLS, identifies shocks recursively, and
% draws the coefficient vector from the plug-in asymptotic distribution
%
%   vec(B_hat) ~ N(vec(B), Sigma_u kron inv(X'X)).
%
% The OLS IRF is the point estimate. Each draw is normalised so the policy
% variable moves by options.desired_policy_impact on impact. The standard
% deviation across draws is the pointwise simulated standard error.
%
% Optional fields:
%   desired_policy_impact   Default: 1
%   draw_covariance_matrix  Default: true
%   confidence_level        Default: 0.95
%
% Requires Statistics and Machine Learning Toolbox for MVNRND and WISHRND.

    if nargin < 6 || isempty(options)
        options = struct();
    end

    if ~isfield(options, 'desired_policy_impact')
        options.desired_policy_impact = 1;
    end

    if ~isfield(options, 'draw_covariance_matrix')
        options.draw_covariance_matrix = true;
    end

    if ~isfield(options, 'confidence_level')
        options.confidence_level = 0.95;
    end

    validateattributes(Y, {'double'}, {'2d', 'nonempty'});
    validateattributes(lags, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(maximum_horizon, {'numeric'}, ...
        {'scalar', 'integer', 'nonnegative'});
    validateattributes(policy_index, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(number_of_draws, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(options.desired_policy_impact, {'numeric'}, ...
        {'scalar', 'finite'});
    validateattributes(options.draw_covariance_matrix, ...
        {'logical', 'numeric'}, {'scalar'});
    validateattributes(options.confidence_level, {'numeric'}, ...
        {'scalar', '>', 0, '<', 1});

    options.draw_covariance_matrix = ...
        logical(options.draw_covariance_matrix);

    %% Validate and trim data

    original_number_of_observations = size(Y, 1);
    number_of_variables = size(Y, 2);

    if policy_index > number_of_variables
        error('policy_index exceeds the number of variables in Y.');
    end

    finite_rows = all(isfinite(Y), 2);
    first_valid_row = find(finite_rows, 1, 'first');
    last_valid_row = find(finite_rows, 1, 'last');

    if isempty(first_valid_row)
        error('Y contains no complete finite observations.');
    end

    Y = Y(first_valid_row:last_valid_row, :);

    if any(~isfinite(Y), 'all')
        error(['Y contains an internal NaN or Inf observation. ', ...
            'Handle internal missing values before estimating the VAR.']);
    end

    number_of_observations = size(Y, 1);
    effective_sample_size = number_of_observations - lags;
    number_of_regressors = 1 + number_of_variables * lags;

    if effective_sample_size <= number_of_regressors
        error(['The effective sample size must exceed the number of ', ...
            'regressors in each VAR equation.']);
    end

    %% Construct VAR regressors

    Y_dependent = Y(lags + 1:end, :);
    X = ones(effective_sample_size, number_of_regressors);

    for current_lag = 1:lags
        columns = 2 + (current_lag - 1) * number_of_variables : ...
            1 + current_lag * number_of_variables;
        X(:, columns) = Y(lags + 1 - current_lag:end - current_lag, :);
    end

    %% OLS estimates

    coefficient_matrix_ols = X \ Y_dependent;
    residuals = Y_dependent - X * coefficient_matrix_ols;

    residual_degrees_of_freedom = ...
        effective_sample_size - number_of_regressors;

    innovation_covariance_ols = ...
        (residuals' * residuals) / residual_degrees_of_freedom;

    innovation_covariance_ols = ...
        make_symmetric_positive_definite(innovation_covariance_ols);

    XX_inverse = (X' * X) \ eye(number_of_regressors);

    coefficient_covariance_asymptotic = kron( ...
        innovation_covariance_ols, XX_inverse);

    coefficient_covariance_asymptotic = ...
        make_symmetric_positive_semidefinite( ...
        coefficient_covariance_asymptotic);

    %% OLS IRF and impact normalisation

    irf_ols_all_shocks = compute_cholesky_irfs( ...
        coefficient_matrix_ols, innovation_covariance_ols, ...
        lags, maximum_horizon);

    irf_ols = squeeze(irf_ols_all_shocks(:, :, policy_index));

    policy_impact_before_normalisation = irf_ols(policy_index, 1);

    if abs(policy_impact_before_normalisation) < 1e-12
        error(['The OLS impact response of the policy variable is ', ...
            'approximately zero, so the IRF cannot be normalised.']);
    end

    normalisation_scaling_factor = options.desired_policy_impact / ...
        policy_impact_before_normalisation;

    irf_ols = irf_ols * normalisation_scaling_factor;

    %% Direct asymptotic coefficient draws using MVNRND

    coefficient_draw_vectors = mvnrnd( ...
        coefficient_matrix_ols(:)', ...
        coefficient_covariance_asymptotic, ...
        number_of_draws);

    irf_draws = NaN(number_of_variables, maximum_horizon + 1, number_of_draws);
    accepted_draws = 0;

    for current_draw = 1:number_of_draws

        coefficient_matrix_draw = reshape( ...
            coefficient_draw_vectors(current_draw, :)', ...
            number_of_regressors, number_of_variables);

        if options.draw_covariance_matrix
            % Plug-in sampling draw centred on Sigma_hat:
            % residual_df * Sigma_hat ~ Wishart(Sigma_u, residual_df).
            innovation_covariance_draw = wishrnd( ...
                innovation_covariance_ols / residual_degrees_of_freedom, ...
                residual_degrees_of_freedom);
            innovation_covariance_draw = ...
                make_symmetric_positive_definite(innovation_covariance_draw);
        else
            innovation_covariance_draw = innovation_covariance_ols;
        end

        irf_draw_all_shocks = compute_cholesky_irfs( ...
            coefficient_matrix_draw, innovation_covariance_draw, ...
            lags, maximum_horizon);

        irf_draw = squeeze(irf_draw_all_shocks(:, :, policy_index));
        policy_impact_draw = irf_draw(policy_index, 1);

        if ~isfinite(policy_impact_draw) || ...
                abs(policy_impact_draw) < 1e-12 || ...
                any(~isfinite(irf_draw), 'all')
            continue
        end

        irf_draw = irf_draw * ...
            (options.desired_policy_impact / policy_impact_draw);

        if any(~isfinite(irf_draw), 'all')
            continue
        end

        accepted_draws = accepted_draws + 1;
        irf_draws(:, :, accepted_draws) = irf_draw;

    end

    irf_draws = irf_draws(:, :, 1:accepted_draws);

    if accepted_draws < max(100, ceil(0.50 * number_of_draws))
        error(['Too few valid Monte Carlo draws were retained. ', ...
            'Inspect the data, lag order and estimated VAR stability.']);
    elseif accepted_draws < number_of_draws
        warning('%d of %d parameter draws produced valid IRFs.', ...
            accepted_draws, number_of_draws);
    end

    %% Pointwise standard errors and intervals

    irf_mean_across_draws = mean(irf_draws, 3);
    irf_pointwise_std = std(irf_draws, 0, 3);

    irf_lower_1se = irf_ols - irf_pointwise_std;
    irf_upper_1se = irf_ols + irf_pointwise_std;

    alpha = 1 - options.confidence_level;
    normal_critical_value = norminv(1 - alpha / 2);

    irf_lower_normal = irf_ols - ...
        normal_critical_value .* irf_pointwise_std;
    irf_upper_normal = irf_ols + ...
        normal_critical_value .* irf_pointwise_std;

    %% Return results

    results = struct();
    results.original_number_of_observations = original_number_of_observations;
    results.first_used_row = first_valid_row;
    results.last_used_row = last_valid_row;
    results.Y = Y;
    results.Y_dependent = Y_dependent;
    results.X = X;
    results.number_of_variables = number_of_variables;
    results.lags = lags;
    results.maximum_horizon = maximum_horizon;
    results.horizons = 0:maximum_horizon;
    results.policy_index = policy_index;
    results.effective_sample_size = effective_sample_size;
    results.number_of_regressors = number_of_regressors;
    results.residual_degrees_of_freedom = residual_degrees_of_freedom;
    results.coefficient_matrix_ols = coefficient_matrix_ols;
    results.residuals = residuals;
    results.innovation_covariance_ols = innovation_covariance_ols;
    results.XX_inverse = XX_inverse;
    results.coefficient_covariance_asymptotic = ...
        coefficient_covariance_asymptotic;
    results.policy_impact_before_normalisation = ...
        policy_impact_before_normalisation;
    results.desired_policy_impact = options.desired_policy_impact;
    results.normalisation_scaling_factor = normalisation_scaling_factor;
    results.irf_ols = irf_ols;
    results.irf_draws = irf_draws;
    results.irf_mean_across_draws = irf_mean_across_draws;
    results.irf_pointwise_std = irf_pointwise_std;
    results.irf_lower_1se = irf_lower_1se;
    results.irf_upper_1se = irf_upper_1se;
    results.confidence_level = options.confidence_level;
    results.normal_critical_value = normal_critical_value;
    results.irf_lower_normal = irf_lower_normal;
    results.irf_upper_normal = irf_upper_normal;
    results.requested_number_of_draws = number_of_draws;
    results.accepted_number_of_draws = accepted_draws;
    results.draw_covariance_matrix = options.draw_covariance_matrix;

end


function irfs = compute_cholesky_irfs( ...
        coefficient_matrix, innovation_covariance, lags, maximum_horizon)
% Output dimensions: response variable x horizon x structural shock.

    number_of_variables = size(innovation_covariance, 1);
    autoregressive_matrices = NaN( ...
        number_of_variables, number_of_variables, lags);

    for current_lag = 1:lags
        rows = 2 + (current_lag - 1) * number_of_variables : ...
            1 + current_lag * number_of_variables;
        autoregressive_matrices(:, :, current_lag) = ...
            coefficient_matrix(rows, :)';
    end

    impact_matrix = chol(innovation_covariance, 'lower');

    moving_average_matrices = zeros( ...
        number_of_variables, number_of_variables, maximum_horizon + 1);
    moving_average_matrices(:, :, 1) = eye(number_of_variables);

    for horizon = 1:maximum_horizon
        current_matrix = zeros(number_of_variables);
        for current_lag = 1:min(lags, horizon)
            current_matrix = current_matrix + ...
                autoregressive_matrices(:, :, current_lag) * ...
                moving_average_matrices(:, :, horizon - current_lag + 1);
        end
        moving_average_matrices(:, :, horizon + 1) = current_matrix;
    end

    irfs = zeros(number_of_variables, maximum_horizon + 1, number_of_variables);

    for horizon = 0:maximum_horizon
        current_irf = moving_average_matrices(:, :, horizon + 1) * impact_matrix;
        irfs(:, horizon + 1, :) = reshape( ...
            current_irf, number_of_variables, 1, number_of_variables);
    end

end


function output_matrix = make_symmetric_positive_definite(input_matrix)

    symmetric_matrix = (input_matrix + input_matrix') / 2;
    [~, flag] = chol(symmetric_matrix);

    if flag == 0
        output_matrix = symmetric_matrix;
        return
    end

    [vectors, values] = eig(symmetric_matrix, 'vector');
    scale = max(1, max(abs(values)));
    values = max(values, 1e-12 * scale);
    output_matrix = vectors * diag(values) * vectors';
    output_matrix = (output_matrix + output_matrix') / 2;

end


function output_matrix = make_symmetric_positive_semidefinite(input_matrix)

    symmetric_matrix = (input_matrix + input_matrix') / 2;
    [vectors, values] = eig(symmetric_matrix, 'vector');
    scale = max(1, max(abs(values)));
    values = max(values, 1e-14 * scale);
    output_matrix = vectors * diag(values) * vectors';
    output_matrix = (output_matrix + output_matrix') / 2;

end
