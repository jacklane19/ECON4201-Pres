function results = var_irf_asymptotic_mvnrnd_seasonal_finite_horizon_v4( ...
        Y, lags, maximum_horizon, policy_index, number_of_draws, options)
%VAR_IRF_ASYMPTOTIC_MVNRND_SEASONAL_FINITE_HORIZON_V4 Frequentist bands for a Cholesky VAR.
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
% Deterministic regressors, including seasonal dummies, are included in
% every VAR equation but are excluded from the dynamic IRF recursion.
%
% Optional fields:
%   desired_policy_impact          Default: 1
%   draw_covariance_matrix         Default: true
%   confidence_level               Default: 0.95
%   deterministic_regressors      Default: no additional regressors
%   deterministic_names           Default: generated variable names
%   response_indices              Default: all response variables
%   maximum_allowed_root          Default: 1.05
%   maximum_allowed_absolute_irf  Default: 100
%   maximum_draw_attempts         Default: max(10000,100*number_of_draws)
%   draw_batch_size               Default: max(100,number_of_draws)
%
% The deterministic_regressors matrix must have the same number of rows as
% Y before leading or trailing incomplete observations are removed.
%
% No stationarity restriction is imposed on the asymptotic parameter draws.
% This is appropriate for a highly persistent levels VAR whose OLS companion
% root may be at or slightly above one. The two editable safeguards reject
% only clearly explosive or numerically extreme finite-horizon draws:
%
%   1. maximum companion root greater than maximum_allowed_root; or
%   2. maximum absolute normalised IRF over the requested horizon greater
%      than maximum_allowed_absolute_irf.
%
% These are loose numerical safeguards, not a claim that retained draws are
% dynamically stable. The routine continues until the requested number of
% valid finite-horizon draws is accepted or maximum_draw_attempts is reached.
%
% Requires Statistics and Machine Learning Toolbox for MVNRND, WISHRND,
% and NORMINV.

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

    if ~isfield(options, 'deterministic_regressors') || ...
            isempty(options.deterministic_regressors)
        options.deterministic_regressors = zeros(size(Y, 1), 0);
    end

    if ~isfield(options, 'deterministic_names')
        options.deterministic_names = strings( ...
            1, size(options.deterministic_regressors, 2));
    end

    if ~isfield(options, 'response_indices') || ...
            isempty(options.response_indices)
        options.response_indices = 1:size(Y, 2);
    end

    if ~isfield(options, 'maximum_allowed_root') || ...
            isempty(options.maximum_allowed_root)
        options.maximum_allowed_root = 1.05;
    end

    if ~isfield(options, 'maximum_allowed_absolute_irf') || ...
            isempty(options.maximum_allowed_absolute_irf)
        options.maximum_allowed_absolute_irf = 100;
    end

    if ~isfield(options, 'maximum_draw_attempts') || ...
            isempty(options.maximum_draw_attempts)
        options.maximum_draw_attempts = max(10000, 100 * number_of_draws);
    end

    if ~isfield(options, 'draw_batch_size') || ...
            isempty(options.draw_batch_size)
        options.draw_batch_size = max(100, number_of_draws);
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
    validateattributes(options.deterministic_regressors, {'numeric'}, ...
        {'2d'});
    validateattributes(options.response_indices, {'numeric'}, ...
        {'vector', 'integer', 'positive'});
    validateattributes(options.maximum_allowed_root, {'numeric'}, ...
        {'scalar', 'finite', '>', 0});
    validateattributes(options.maximum_allowed_absolute_irf, {'numeric'}, ...
        {'scalar', 'finite', '>', 0});
    validateattributes(options.maximum_draw_attempts, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(options.draw_batch_size, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});

    if options.maximum_draw_attempts < number_of_draws
        error(['maximum_draw_attempts must be at least as large as ', ...
            'number_of_draws.']);
    end

    options.draw_covariance_matrix = ...
        logical(options.draw_covariance_matrix);
    options.response_indices = options.response_indices(:)';

    %% Validate and align data

    original_number_of_observations = size(Y, 1);
    number_of_variables = size(Y, 2);

    if policy_index > number_of_variables
        error('policy_index exceeds the number of variables in Y.');
    end

    if any(options.response_indices > number_of_variables) || ...
            numel(unique(options.response_indices)) ~= ...
            numel(options.response_indices)
        error(['response_indices must contain unique valid indices ', ...
            'for variables in Y.']);
    end

    deterministic_regressors = options.deterministic_regressors;

    if size(deterministic_regressors, 1) ~= ...
            original_number_of_observations
        error(['deterministic_regressors must have the same number ', ...
            'of rows as Y.']);
    end

    if any(~isfinite(deterministic_regressors), 'all')
        error('deterministic_regressors contains NaN or Inf values.');
    end

    number_of_deterministic_regressors = ...
        size(deterministic_regressors, 2);

    deterministic_names = string(options.deterministic_names);
    deterministic_names = deterministic_names(:)';

    if isempty(deterministic_names) && ...
            number_of_deterministic_regressors > 0
        deterministic_names = "Deterministic" + ...
            string(1:number_of_deterministic_regressors);
    elseif numel(deterministic_names) ~= ...
            number_of_deterministic_regressors
        error(['deterministic_names must have one entry for each ', ...
            'column of deterministic_regressors.']);
    end

    finite_rows = all(isfinite(Y), 2);
    first_valid_row = find(finite_rows, 1, 'first');
    last_valid_row = find(finite_rows, 1, 'last');

    if isempty(first_valid_row)
        error('Y contains no complete finite observations.');
    end

    Y = Y(first_valid_row:last_valid_row, :);
    deterministic_regressors = deterministic_regressors( ...
        first_valid_row:last_valid_row, :);

    if any(~isfinite(Y), 'all')
        error(['Y contains an internal NaN or Inf observation. ', ...
            'Handle internal missing values before estimating the VAR.']);
    end

    number_of_observations = size(Y, 1);
    effective_sample_size = number_of_observations - lags;
    number_of_regressors = 1 + ...
        number_of_deterministic_regressors + number_of_variables * lags;

    if effective_sample_size <= number_of_regressors
        error(['The effective sample size must exceed the number of ', ...
            'regressors in each VAR equation.']);
    end

    %% Construct VAR regressors

    Y_dependent = Y(lags + 1:end, :);
    deterministic_estimation_sample = ...
        deterministic_regressors(lags + 1:end, :);

    X = [ ...
        ones(effective_sample_size, 1), ...
        deterministic_estimation_sample, ...
        zeros(effective_sample_size, number_of_variables * lags)];

    first_lag_column = 2 + number_of_deterministic_regressors;

    for current_lag = 1:lags
        columns = first_lag_column + ...
            (current_lag - 1) * number_of_variables : ...
            first_lag_column - 1 + ...
            current_lag * number_of_variables;

        X(:, columns) = Y( ...
            lags + 1 - current_lag:end - current_lag, :);
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

    %% OLS companion-root diagnostic

    maximum_companion_root_ols = compute_maximum_companion_root( ...
        coefficient_matrix_ols, lags, number_of_variables, ...
        number_of_deterministic_regressors);

    ols_is_dynamically_stable = ...
        isfinite(maximum_companion_root_ols) && ...
        maximum_companion_root_ols < 1;

    ols_exceeds_numerical_root_safeguard = ...
        ~isfinite(maximum_companion_root_ols) || ...
        maximum_companion_root_ols > options.maximum_allowed_root;

    if maximum_companion_root_ols >= 1 && ...
            ~ols_exceeds_numerical_root_safeguard
        fprintf(['Note: the estimated OLS VAR has a largest absolute ', ...
            'companion root of %.6f. No stationarity filter is imposed ', ...
            'on the asymptotic draws because this is a persistent levels ', ...
            'VAR. Only roots above the loose numerical safeguard of ', ...
            '%.6f, or numerically extreme finite-horizon IRFs, are ', ...
            'rejected.\n'], ...
            maximum_companion_root_ols, options.maximum_allowed_root);
    end

    if ols_exceeds_numerical_root_safeguard
        warning(['The OLS companion root itself exceeds the loose ', ...
            'numerical root safeguard. The OLS IRF remains the point ', ...
            'estimate, but the simulation distribution is truncated at ', ...
            'a maximum root of %.6f. Consider increasing ', ...
            'maximum_allowed_root after inspecting the finite-horizon ', ...
            'responses.'], options.maximum_allowed_root);
    end

    %% OLS IRF and impact normalisation

    irf_ols_all_shocks = compute_cholesky_irfs( ...
        coefficient_matrix_ols, innovation_covariance_ols, ...
        lags, maximum_horizon, number_of_deterministic_regressors);

    irf_ols_full = squeeze( ...
        irf_ols_all_shocks(:, :, policy_index));

    policy_impact_before_normalisation = ...
        irf_ols_full(policy_index, 1);

    if abs(policy_impact_before_normalisation) < 1e-12
        error(['The OLS impact response of the policy variable is ', ...
            'approximately zero, so the IRF cannot be normalised.']);
    end

    normalisation_scaling_factor = options.desired_policy_impact / ...
        policy_impact_before_normalisation;

    irf_ols_full = irf_ols_full * normalisation_scaling_factor;
    irf_ols = irf_ols_full(options.response_indices, :);

    %% Direct asymptotic coefficient draws using MVNRND
    %
    % No strict stationarity filter is applied. Draws are rejected only if
    % they are nonfinite, have a companion root above the loose numerical
    % safeguard, cannot be normalised, or produce an extreme normalised IRF
    % over the requested finite horizon.

    number_of_returned_responses = numel(options.response_indices);
    irf_draws = NaN( ...
        number_of_returned_responses, ...
        maximum_horizon + 1, ...
        number_of_draws);

    accepted_draws = 0;
    attempted_draws = 0;
    rejected_root_safeguard_draws = 0;
    rejected_extreme_irf_draws = 0;
    rejected_invalid_draws = 0;
    maximum_roots_accepted = NaN(number_of_draws, 1);
    maximum_absolute_irfs_accepted = NaN(number_of_draws, 1);

    while accepted_draws < number_of_draws && ...
            attempted_draws < options.maximum_draw_attempts

        available_attempts = ...
            options.maximum_draw_attempts - attempted_draws;
        current_batch_size = min(options.draw_batch_size, available_attempts);

        coefficient_draw_vectors = mvnrnd( ...
            coefficient_matrix_ols(:)', ...
            coefficient_covariance_asymptotic, ...
            current_batch_size);

        for draw_in_batch = 1:current_batch_size

            attempted_draws = attempted_draws + 1;

            coefficient_matrix_draw = reshape( ...
                coefficient_draw_vectors(draw_in_batch, :)', ...
                number_of_regressors, number_of_variables);

            maximum_companion_root_draw = compute_maximum_companion_root( ...
                coefficient_matrix_draw, lags, number_of_variables, ...
                number_of_deterministic_regressors);

            if ~isfinite(maximum_companion_root_draw)
                rejected_invalid_draws = rejected_invalid_draws + 1;
                continue
            end

            if maximum_companion_root_draw > options.maximum_allowed_root
                rejected_root_safeguard_draws = ...
                    rejected_root_safeguard_draws + 1;
                continue
            end

            if options.draw_covariance_matrix
                % Plug-in sampling draw centred on Sigma_hat:
                % E[Sigma_draw | Sigma_hat] = Sigma_hat.
                innovation_covariance_draw = wishrnd( ...
                    innovation_covariance_ols / ...
                    residual_degrees_of_freedom, ...
                    residual_degrees_of_freedom);
                innovation_covariance_draw = ...
                    make_symmetric_positive_definite( ...
                    innovation_covariance_draw);
            else
                innovation_covariance_draw = innovation_covariance_ols;
            end

            irf_draw_all_shocks = compute_cholesky_irfs( ...
                coefficient_matrix_draw, innovation_covariance_draw, ...
                lags, maximum_horizon, ...
                number_of_deterministic_regressors);

            irf_draw_full = squeeze( ...
                irf_draw_all_shocks(:, :, policy_index));
            policy_impact_draw = irf_draw_full(policy_index, 1);

            if ~isfinite(policy_impact_draw) || ...
                    abs(policy_impact_draw) < 1e-12 || ...
                    any(~isfinite(irf_draw_full), 'all')
                rejected_invalid_draws = rejected_invalid_draws + 1;
                continue
            end

            irf_draw_full = irf_draw_full * ...
                (options.desired_policy_impact / policy_impact_draw);

            if any(~isfinite(irf_draw_full), 'all')
                rejected_invalid_draws = rejected_invalid_draws + 1;
                continue
            end

            maximum_absolute_irf_draw = ...
                max(abs(irf_draw_full), [], 'all');

            if maximum_absolute_irf_draw > ...
                    options.maximum_allowed_absolute_irf
                rejected_extreme_irf_draws = ...
                    rejected_extreme_irf_draws + 1;
                continue
            end

            accepted_draws = accepted_draws + 1;
            irf_draws(:, :, accepted_draws) = ...
                irf_draw_full(options.response_indices, :);
            maximum_roots_accepted(accepted_draws) = ...
                maximum_companion_root_draw;
            maximum_absolute_irfs_accepted(accepted_draws) = ...
                maximum_absolute_irf_draw;

            if accepted_draws == number_of_draws
                break
            end

        end
    end

    if accepted_draws < number_of_draws
        error(['Only %d valid finite-horizon draws were accepted after ', ...
            '%d attempts. %d draws exceeded the loose root safeguard ', ...
            '(root > %.6f), %d exceeded the absolute-IRF safeguard ', ...
            '(maximum absolute IRF > %.6f), and %d were invalid. ', ...
            'Increase maximum_draw_attempts or inspect and cautiously ', ...
            'relax the numerical safeguards.'], ...
            accepted_draws, attempted_draws, ...
            rejected_root_safeguard_draws, ...
            options.maximum_allowed_root, ...
            rejected_extreme_irf_draws, ...
            options.maximum_allowed_absolute_irf, ...
            rejected_invalid_draws);
    end

    maximum_roots_accepted = ...
        maximum_roots_accepted(1:accepted_draws);
    maximum_absolute_irfs_accepted = ...
        maximum_absolute_irfs_accepted(1:accepted_draws);
    draw_acceptance_rate = accepted_draws / attempted_draws;

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
    results.original_number_of_observations = ...
        original_number_of_observations;
    results.first_used_row = first_valid_row;
    results.last_used_row = last_valid_row;
    results.Y = Y;
    results.Y_dependent = Y_dependent;
    results.X = X;
    results.deterministic_regressors = deterministic_regressors;
    results.deterministic_estimation_sample = ...
        deterministic_estimation_sample;
    results.deterministic_names = deterministic_names;
    results.number_of_deterministic_regressors = ...
        number_of_deterministic_regressors;
    results.number_of_variables = number_of_variables;
    results.response_indices = options.response_indices;
    results.lags = lags;
    results.maximum_horizon = maximum_horizon;
    results.horizons = 0:maximum_horizon;
    results.policy_index = policy_index;
    results.effective_sample_size = effective_sample_size;
    results.number_of_regressors = number_of_regressors;
    results.residual_degrees_of_freedom = ...
        residual_degrees_of_freedom;
    results.coefficient_matrix_ols = coefficient_matrix_ols;
    results.residuals = residuals;
    results.innovation_covariance_ols = innovation_covariance_ols;
    results.XX_inverse = XX_inverse;
    results.coefficient_covariance_asymptotic = ...
        coefficient_covariance_asymptotic;
    results.policy_impact_before_normalisation = ...
        policy_impact_before_normalisation;
    results.desired_policy_impact = options.desired_policy_impact;
    results.normalisation_scaling_factor = ...
        normalisation_scaling_factor;
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
    results.attempted_number_of_draws = attempted_draws;
    results.rejected_root_safeguard_draws = ...
        rejected_root_safeguard_draws;
    results.rejected_extreme_irf_draws = rejected_extreme_irf_draws;
    results.rejected_invalid_draws = rejected_invalid_draws;
    results.draw_acceptance_rate = draw_acceptance_rate;
    results.maximum_allowed_root = options.maximum_allowed_root;
    results.maximum_allowed_absolute_irf = ...
        options.maximum_allowed_absolute_irf;
    results.maximum_draw_attempts = options.maximum_draw_attempts;
    results.draw_batch_size = options.draw_batch_size;
    results.maximum_companion_root_ols = maximum_companion_root_ols;
    results.ols_is_dynamically_stable = ols_is_dynamically_stable;
    results.ols_exceeds_numerical_root_safeguard = ...
        ols_exceeds_numerical_root_safeguard;
    results.maximum_companion_roots_accepted = maximum_roots_accepted;
    results.maximum_absolute_irfs_accepted = ...
        maximum_absolute_irfs_accepted;
    results.draw_covariance_matrix = ...
        options.draw_covariance_matrix;

end


function irfs = compute_cholesky_irfs( ...
        coefficient_matrix, innovation_covariance, lags, ...
        maximum_horizon, number_of_deterministic_regressors)
% Output dimensions: response variable x horizon x structural shock.
% Constants and deterministic regressors are excluded from the recursion.

    number_of_variables = size(innovation_covariance, 1);
    autoregressive_matrices = NaN( ...
        number_of_variables, number_of_variables, lags);

    first_lag_row = 2 + number_of_deterministic_regressors;

    for current_lag = 1:lags
        rows = first_lag_row + ...
            (current_lag - 1) * number_of_variables : ...
            first_lag_row - 1 + current_lag * number_of_variables;

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
                moving_average_matrices( ...
                :, :, horizon - current_lag + 1);
        end

        moving_average_matrices(:, :, horizon + 1) = current_matrix;
    end

    irfs = zeros( ...
        number_of_variables, maximum_horizon + 1, number_of_variables);

    for horizon = 0:maximum_horizon
        current_irf = ...
            moving_average_matrices(:, :, horizon + 1) * impact_matrix;
        irfs(:, horizon + 1, :) = reshape( ...
            current_irf, number_of_variables, 1, number_of_variables);
    end

end


function maximum_root = compute_maximum_companion_root( ...
        coefficient_matrix, lags, number_of_variables, ...
        number_of_deterministic_regressors)
%COMPUTE_MAXIMUM_COMPANION_ROOT Largest absolute VAR companion eigenvalue.
% Constants and deterministic regressors are excluded from the companion
% matrix because they do not affect dynamic stability.

    first_lag_row = 2 + number_of_deterministic_regressors;
    top_companion_block = zeros( ...
        number_of_variables, number_of_variables * lags);

    for current_lag = 1:lags
        rows = first_lag_row + ...
            (current_lag - 1) * number_of_variables : ...
            first_lag_row - 1 + current_lag * number_of_variables;

        columns = (current_lag - 1) * number_of_variables + 1 : ...
            current_lag * number_of_variables;

        top_companion_block(:, columns) = ...
            coefficient_matrix(rows, :)';
    end

    if lags == 1
        companion_matrix = top_companion_block;
    else
        companion_matrix = [ ...
            top_companion_block; ...
            eye(number_of_variables * (lags - 1)), ...
            zeros(number_of_variables * (lags - 1), ...
            number_of_variables)];
    end

    companion_eigenvalues = eig(companion_matrix);
    maximum_root = max(abs(companion_eigenvalues));

end


function output_matrix = ...
        make_symmetric_positive_definite(input_matrix)

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


function output_matrix = ...
        make_symmetric_positive_semidefinite(input_matrix)

    symmetric_matrix = (input_matrix + input_matrix') / 2;
    [vectors, values] = eig(symmetric_matrix, 'vector');
    scale = max(1, max(abs(values)));
    values = max(values, 0);

    % Remove tiny numerical negative eigenvalues and add minimal jitter only
    % when needed by MVNRND's covariance-factorisation step.
    output_matrix = vectors * diag(values) * vectors';
    output_matrix = (output_matrix + output_matrix') / 2;

    [~, flag] = chol(output_matrix + 1e-14 * scale * eye(size(output_matrix)));
    if flag ~= 0
        output_matrix = output_matrix + ...
            1e-12 * scale * eye(size(output_matrix));
    end

end
