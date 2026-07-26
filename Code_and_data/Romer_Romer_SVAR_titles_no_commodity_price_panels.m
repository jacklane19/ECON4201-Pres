%% Four Cholesky SVAR specifications with model-specific variable tags
%
% Model suffixes:
%   _mpshock      = Romer-Romer monetary policy shock, no commodity prices
%   _mpshock_com  = Romer-Romer monetary policy shock, with commodity prices
%   _ff           = Federal funds rate, no commodity prices
%   _ff_com       = Federal funds rate, with commodity prices
%
% Every model-specific object uses its model suffix so that estimates,
% impulse responses, normalisations, horizons, figures and axes are retained
% separately in the MATLAB workspace.

clear;
clc;
close all;

%% Add the Empirical Macro Toolbox to the MATLAB path

addpath(fullfile(pwd, '..', 'BVAR_', 'bvartools'));

%% Import the Romer and Romer data

data_table = readtable( ...
    "RomerandRomerDataAppendix.xls", ...
    'Sheet', "DATA BY MONTH");

% Display imported variable names
disp(data_table.Properties.VariableNames);

%% Extract the common data series

% Natural log of industrial production
industrial_production_log = data_table.LNIPNSA;

% Natural log of producer prices
producer_price_log = data_table.LNPPINSA;

% Natural log of world commodity prices
commodity_prices_log = data_table.LNWCP;

% Cumulated Romer and Romer monetary policy shock
monetary_shock = data_table.SUMSHCK;

% Federal funds rate, measured in percentage points
fed_funds = data_table.FF;

%% Check whether SUMSHCK and SUMSHCKF are identical

same_shock_series = isequal( ...
    data_table.SUMSHCK, ...
    data_table.SUMSHCKF);

fprintf('Are SUMSHCK and SUMSHCKF identical? %d\n', ...
    same_shock_series);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Monetary policy shock VAR without commodity prices: _mpshock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct the VAR dataset

% Cholesky ordering:
% 1. Log industrial production
% 2. Log producer price level
% 3. Cumulated Romer-Romer monetary policy shock

y_mpshock = [
    industrial_production_log, ...
    producer_price_log, ...
    monetary_shock
    ];

%% Estimate the Cholesky SVAR

lags_mpshock = 36;

opt_mpshock.hor = 48;
opt_mpshock.K   = 0;

VAR_mpshock = cvar_( ...
    y_mpshock, ...
    lags_mpshock, ...
    opt_mpshock);

%% Extract responses to the monetary policy shock

monetary_shock_index_mpshock = 3;
number_of_variables_mpshock = 3;

irfs_mpshock = squeeze( ...
    VAR_mpshock.ir_ols(:, :, monetary_shock_index_mpshock, :) ...
    );

% Ensure rows represent variables and columns represent horizons
if size(irfs_mpshock, 1) ~= number_of_variables_mpshock && ...
        size(irfs_mpshock, 2) == number_of_variables_mpshock
    irfs_mpshock = irfs_mpshock';
end

if size(irfs_mpshock, 1) ~= number_of_variables_mpshock
    error(['The extracted IRF matrix for the _mpshock model does not ', ...
        'have three variables in its rows. Check VAR_mpshock.ir_ols.']);
end

%% Normalise to a one-percentage-point monetary policy innovation

policy_variable_index_mpshock = 3;

policy_impact_response_mpshock = ...
    irfs_mpshock(policy_variable_index_mpshock, 1);

if abs(policy_impact_response_mpshock) < 1e-12
    error(['The impact response of SUMSHCK in the _mpshock model is ', ...
        'approximately zero, so the IRFs cannot be normalised.']);
end

desired_policy_impact_mpshock = 1;

scaling_factor_mpshock = ...
    desired_policy_impact_mpshock / policy_impact_response_mpshock;

irfs_rescaled_mpshock = ...
    irfs_mpshock * scaling_factor_mpshock;

%% Convert logged responses into percentages

irfs_for_plot_mpshock = irfs_rescaled_mpshock;

% Industrial production and producer prices are in natural logarithms.
irfs_for_plot_mpshock(1:2, :) = ...
    100 * irfs_rescaled_mpshock(1:2, :);

% SUMSHCK remains in percentage points.
irfs_for_plot_mpshock(3, :) = ...
    irfs_rescaled_mpshock(3, :);

%% Display normalisation information

fprintf('\n_mpshock original policy impact response: %.6f\n', ...
    policy_impact_response_mpshock);

fprintf('_mpshock scaling factor: %.6f\n', ...
    scaling_factor_mpshock);

fprintf(['_mpshock rescaled policy impact response ', ...
    '(should equal 1): %.6f\n'], ...
    irfs_for_plot_mpshock(policy_variable_index_mpshock, 1));

%% Construct the horizon vector

number_of_horizons_mpshock = ...
    size(irfs_for_plot_mpshock, 2);

horizons_mpshock = ...
    0:(number_of_horizons_mpshock - 1);

%% Plot the _mpshock impulse responses

panel_indices_mpshock = [3, 1, 2];

panel_titles_mpshock = {
    'Effect on the Cumulated Shock'
    'Effect on Output'
    'Effect on the Price Level'
    };

y_axis_labels_mpshock = {
    'Percentage Points'
    'Percent'
    'Percent'
    };

y_axis_limits_mpshock = {
    [0, 1.2]
    [-5, 2]
    [-7, 1]
    };

y_axis_ticks_mpshock = {
    0:0.2:1.2
    -5:1:2
    -7:1:1
    };

figure_title_mpshock = ...
    'Romer-Romer Monetary Policy Shock VAR';

[figure_mpshock, plot_layout_mpshock, axes_mpshock] = ...
    plot_tagged_irfs( ...
    irfs_for_plot_mpshock, ...
    horizons_mpshock, ...
    panel_indices_mpshock, ...
    panel_titles_mpshock, ...
    y_axis_labels_mpshock, ...
    y_axis_limits_mpshock, ...
    y_axis_ticks_mpshock, ...
    [100, 50, 750, 1050], ...
    figure_title_mpshock);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Monetary policy shock VAR with commodity prices: _mpshock_com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct the VAR dataset

% Cholesky ordering:
% 1. Log industrial production
% 2. Log producer price level
% 3. Log world commodity prices
% 4. Cumulated Romer-Romer monetary policy shock

y_mpshock_com = [
    industrial_production_log, ...
    producer_price_log, ...
    commodity_prices_log, ...
    monetary_shock
    ];

%% Estimate the Cholesky SVAR

lags_mpshock_com = 36;

opt_mpshock_com.hor = 48;
opt_mpshock_com.K   = 0;

VAR_mpshock_com = cvar_( ...
    y_mpshock_com, ...
    lags_mpshock_com, ...
    opt_mpshock_com);

%% Extract responses to the monetary policy shock

monetary_shock_index_mpshock_com = 4;
number_of_variables_mpshock_com = 4;

irfs_mpshock_com = squeeze( ...
    VAR_mpshock_com.ir_ols(:, :, monetary_shock_index_mpshock_com, :) ...
    );

% Ensure rows represent variables and columns represent horizons
if size(irfs_mpshock_com, 1) ~= number_of_variables_mpshock_com && ...
        size(irfs_mpshock_com, 2) == number_of_variables_mpshock_com
    irfs_mpshock_com = irfs_mpshock_com';
end

if size(irfs_mpshock_com, 1) ~= number_of_variables_mpshock_com
    error(['The extracted IRF matrix for the _mpshock_com model does ', ...
        'not have four variables in its rows. ', ...
        'Check VAR_mpshock_com.ir_ols.']);
end

%% Normalise to a one-percentage-point monetary policy innovation

policy_variable_index_mpshock_com = 4;

policy_impact_response_mpshock_com = ...
    irfs_mpshock_com(policy_variable_index_mpshock_com, 1);

if abs(policy_impact_response_mpshock_com) < 1e-12
    error(['The impact response of SUMSHCK in the _mpshock_com model ', ...
        'is approximately zero, so the IRFs cannot be normalised.']);
end

desired_policy_impact_mpshock_com = 1;

scaling_factor_mpshock_com = ...
    desired_policy_impact_mpshock_com / ...
    policy_impact_response_mpshock_com;

irfs_rescaled_mpshock_com = ...
    irfs_mpshock_com * scaling_factor_mpshock_com;

%% Convert logged responses into percentages

irfs_for_plot_mpshock_com = irfs_rescaled_mpshock_com;

% Industrial production, producer prices and commodity prices are logs.
irfs_for_plot_mpshock_com(1:3, :) = ...
    100 * irfs_rescaled_mpshock_com(1:3, :);

% SUMSHCK remains in percentage points.
irfs_for_plot_mpshock_com(4, :) = ...
    irfs_rescaled_mpshock_com(4, :);

%% Display normalisation information

fprintf('\n_mpshock_com original policy impact response: %.6f\n', ...
    policy_impact_response_mpshock_com);

fprintf('_mpshock_com scaling factor: %.6f\n', ...
    scaling_factor_mpshock_com);

fprintf(['_mpshock_com rescaled policy impact response ', ...
    '(should equal 1): %.6f\n'], ...
    irfs_for_plot_mpshock_com(policy_variable_index_mpshock_com, 1));

%% Construct the horizon vector

number_of_horizons_mpshock_com = ...
    size(irfs_for_plot_mpshock_com, 2);

horizons_mpshock_com = ...
    0:(number_of_horizons_mpshock_com - 1);

%% Plot the _mpshock_com impulse responses

% Commodity prices remain in the estimated VAR but their IRF is not plotted.
panel_indices_mpshock_com = [4, 1, 2];

panel_titles_mpshock_com = {
    'Effect on the Cumulated Shock'
    'Effect on Output'
    'Effect on the Price Level'
    };

y_axis_labels_mpshock_com = {
    'Percentage Points'
    'Percent'
    'Percent'
    };

y_axis_limits_mpshock_com = {
    [0, 1.2]
    [-5, 2]
    [-7, 1]
    };

y_axis_ticks_mpshock_com = {
    0:0.2:1.2
    -5:1:2
    -7:1:1
    };

figure_title_mpshock_com = ...
    'Romer-Romer Monetary Policy Shock VAR with Commodity Prices';

[figure_mpshock_com, plot_layout_mpshock_com, axes_mpshock_com] = ...
    plot_tagged_irfs( ...
    irfs_for_plot_mpshock_com, ...
    horizons_mpshock_com, ...
    panel_indices_mpshock_com, ...
    panel_titles_mpshock_com, ...
    y_axis_labels_mpshock_com, ...
    y_axis_limits_mpshock_com, ...
    y_axis_ticks_mpshock_com, ...
    [900, 50, 800, 1050], ...
    figure_title_mpshock_com);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Federal funds rate VAR without commodity prices: _ff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct the VAR dataset

% Cholesky ordering:
% 1. Log industrial production
% 2. Log producer price level
% 3. Federal funds rate

y_ff = [
    industrial_production_log, ...
    producer_price_log, ...
    fed_funds
    ];

%% Estimate the Cholesky SVAR

lags_ff = 36;

opt_ff.hor = 48;
opt_ff.K   = 0;

VAR_ff = cvar_( ...
    y_ff, ...
    lags_ff, ...
    opt_ff);

%% Extract responses to the federal funds rate shock

monetary_shock_index_ff = 3;
number_of_variables_ff = 3;

irfs_ff = squeeze( ...
    VAR_ff.ir_ols(:, :, monetary_shock_index_ff, :) ...
    );

% Ensure rows represent variables and columns represent horizons
if size(irfs_ff, 1) ~= number_of_variables_ff && ...
        size(irfs_ff, 2) == number_of_variables_ff
    irfs_ff = irfs_ff';
end

if size(irfs_ff, 1) ~= number_of_variables_ff
    error(['The extracted IRF matrix for the _ff model does not have ', ...
        'three variables in its rows. Check VAR_ff.ir_ols.']);
end

%% Normalise to a one-percentage-point federal funds rate innovation

policy_variable_index_ff = 3;

policy_impact_response_ff = ...
    irfs_ff(policy_variable_index_ff, 1);

if abs(policy_impact_response_ff) < 1e-12
    error(['The impact response of the federal funds rate in the _ff ', ...
        'model is approximately zero, so the IRFs cannot be normalised.']);
end

desired_policy_impact_ff = 1;

scaling_factor_ff = ...
    desired_policy_impact_ff / policy_impact_response_ff;

irfs_rescaled_ff = ...
    irfs_ff * scaling_factor_ff;

%% Convert logged responses into percentages

irfs_for_plot_ff = irfs_rescaled_ff;

% Industrial production and producer prices are in natural logarithms.
irfs_for_plot_ff(1:2, :) = ...
    100 * irfs_rescaled_ff(1:2, :);

% The federal funds rate remains in percentage points.
irfs_for_plot_ff(3, :) = ...
    irfs_rescaled_ff(3, :);

%% Display normalisation information

fprintf('\n_ff original policy impact response: %.6f\n', ...
    policy_impact_response_ff);

fprintf('_ff scaling factor: %.6f\n', ...
    scaling_factor_ff);

fprintf(['_ff rescaled policy impact response ', ...
    '(should equal 1): %.6f\n'], ...
    irfs_for_plot_ff(policy_variable_index_ff, 1));

%% Construct the horizon vector

number_of_horizons_ff = ...
    size(irfs_for_plot_ff, 2);

horizons_ff = ...
    0:(number_of_horizons_ff - 1);

%% Plot the _ff impulse responses

panel_indices_ff = [3, 1, 2];

panel_titles_ff = {
    'Effect on the Federal Funds Rate'
    'Effect on Output'
    'Effect on the Price Level'
    };

y_axis_labels_ff = {
    'Percentage Points'
    'Percent'
    'Percent'
    };

y_axis_limits_ff = {
    [0, 1.2]
    [-5, 2]
    [-7, 1]
    };

y_axis_ticks_ff = {
    0:0.2:1.2
    -5:1:2
    -7:1:1
    };

figure_title_ff = ...
    'Federal Funds Rate VAR';

[figure_ff, plot_layout_ff, axes_ff] = ...
    plot_tagged_irfs( ...
    irfs_for_plot_ff, ...
    horizons_ff, ...
    panel_indices_ff, ...
    panel_titles_ff, ...
    y_axis_labels_ff, ...
    y_axis_limits_ff, ...
    y_axis_ticks_ff, ...
    [100, 100, 750, 1050], ...
    figure_title_ff);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Federal funds rate VAR with commodity prices: _ff_com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct the VAR dataset

% Cholesky ordering:
% 1. Log industrial production
% 2. Log producer price level
% 3. Log world commodity prices
% 4. Federal funds rate

y_ff_com = [
    industrial_production_log, ...
    producer_price_log, ...
    commodity_prices_log, ...
    fed_funds
    ];

%% Estimate the Cholesky SVAR

lags_ff_com = 36;

opt_ff_com.hor = 48;
opt_ff_com.K   = 0;

VAR_ff_com = cvar_( ...
    y_ff_com, ...
    lags_ff_com, ...
    opt_ff_com);

%% Extract responses to the federal funds rate shock

monetary_shock_index_ff_com = 4;
number_of_variables_ff_com = 4;

irfs_ff_com = squeeze( ...
    VAR_ff_com.ir_ols(:, :, monetary_shock_index_ff_com, :) ...
    );

% Ensure rows represent variables and columns represent horizons
if size(irfs_ff_com, 1) ~= number_of_variables_ff_com && ...
        size(irfs_ff_com, 2) == number_of_variables_ff_com
    irfs_ff_com = irfs_ff_com';
end

if size(irfs_ff_com, 1) ~= number_of_variables_ff_com
    error(['The extracted IRF matrix for the _ff_com model does not ', ...
        'have four variables in its rows. Check VAR_ff_com.ir_ols.']);
end

%% Normalise to a one-percentage-point federal funds rate innovation

policy_variable_index_ff_com = 4;

policy_impact_response_ff_com = ...
    irfs_ff_com(policy_variable_index_ff_com, 1);

if abs(policy_impact_response_ff_com) < 1e-12
    error(['The impact response of the federal funds rate in the ', ...
        '_ff_com model is approximately zero, so the IRFs cannot ', ...
        'be normalised.']);
end

desired_policy_impact_ff_com = 1;

scaling_factor_ff_com = ...
    desired_policy_impact_ff_com / policy_impact_response_ff_com;

irfs_rescaled_ff_com = ...
    irfs_ff_com * scaling_factor_ff_com;

%% Convert logged responses into percentages

irfs_for_plot_ff_com = irfs_rescaled_ff_com;

% Industrial production, producer prices and commodity prices are logs.
irfs_for_plot_ff_com(1:3, :) = ...
    100 * irfs_rescaled_ff_com(1:3, :);

% The federal funds rate remains in percentage points.
irfs_for_plot_ff_com(4, :) = ...
    irfs_rescaled_ff_com(4, :);

%% Display normalisation information

fprintf('\n_ff_com original policy impact response: %.6f\n', ...
    policy_impact_response_ff_com);

fprintf('_ff_com scaling factor: %.6f\n', ...
    scaling_factor_ff_com);

fprintf(['_ff_com rescaled policy impact response ', ...
    '(should equal 1): %.6f\n'], ...
    irfs_for_plot_ff_com(policy_variable_index_ff_com, 1));

%% Construct the horizon vector

number_of_horizons_ff_com = ...
    size(irfs_for_plot_ff_com, 2);

horizons_ff_com = ...
    0:(number_of_horizons_ff_com - 1);

%% Plot the _ff_com impulse responses

% Commodity prices remain in the estimated VAR but their IRF is not plotted.
panel_indices_ff_com = [4, 1, 2];

panel_titles_ff_com = {
    'Effect on the Federal Funds Rate'
    'Effect on Output'
    'Effect on the Price Level'
    };

y_axis_labels_ff_com = {
    'Percentage Points'
    'Percent'
    'Percent'
    };

y_axis_limits_ff_com = {
    [0, 1.2]
    [-5, 2]
    [-7, 1]
    };

y_axis_ticks_ff_com = {
    0:0.2:1.2
    -5:1:2
    -7:1:1
    };

figure_title_ff_com = ...
    'Federal Funds Rate VAR with Commodity Prices';

[figure_ff_com, plot_layout_ff_com, axes_ff_com] = ...
    plot_tagged_irfs( ...
    irfs_for_plot_ff_com, ...
    horizons_ff_com, ...
    panel_indices_ff_com, ...
    panel_titles_ff_com, ...
    y_axis_labels_ff_com, ...
    y_axis_limits_ff_com, ...
    y_axis_ticks_ff_com, ...
    [900, 100, 800, 1050], ...
    figure_title_ff_com);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local plotting function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [figure_handle, layout_handle, axes_handles] = plot_tagged_irfs( ...
        irfs_for_plot, ...
        horizons, ...
        panel_indices, ...
        panel_titles, ...
        y_axis_labels, ...
        y_axis_limits, ...
        y_axis_ticks, ...
        figure_position, ...
        figure_title)
%PLOT_TAGGED_IRFS Plot one specification's impulse response functions.
%
% This function uses a local workspace, so its temporary plotting variables
% cannot overwrite any of the tagged variables in the main script workspace.

    number_of_panels = numel(panel_indices);
    panel_letters = 'abcdefghijklmnopqrstuvwxyz';

    figure_handle = figure( ...
        'Position', figure_position, ...
        'Color', 'w');

    layout_handle = tiledlayout(number_of_panels, 1);
    layout_handle.TileSpacing = 'compact';
    layout_handle.Padding = 'compact';

    % Add the model-specific title above the full tiled figure.
    title( ...
        layout_handle, ...
        figure_title, ...
        'FontSize', 20, ...
        'FontWeight', 'bold');

    axes_handles = gobjects(number_of_panels, 1);

    for panel = 1:number_of_panels

        axes_handles(panel) = nexttile(layout_handle);
        current_axis = axes_handles(panel);
        current_variable = panel_indices(panel);

        plot( ...
            current_axis, ...
            horizons, ...
            irfs_for_plot(current_variable, :), ...
            'k-', ...
            'LineWidth', 1.5);

        hold(current_axis, 'on');

        yline( ...
            current_axis, ...
            0, ...
            'k:', ...
            'LineWidth', 1);

        title( ...
            current_axis, ...
            sprintf('%s. %s', ...
            panel_letters(panel), ...
            panel_titles{panel}), ...
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

        xlim(current_axis, [0, 48]);
        xticks(current_axis, 0:3:48);

        current_axis.FontSize = 14;
        current_axis.Box = 'on';
        current_axis.LineWidth = 1;
    end

    linkaxes(axes_handles, 'x');
end


