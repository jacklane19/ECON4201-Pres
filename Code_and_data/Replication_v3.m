%% Cholesky SVAR using the Romer and Romer monetary policy shock
% IRFs are normalised to a one-percentage-point innovation in SUMSHCK.
% Logged output and price responses are converted into percentages.
% The plots use the same y-axis limits as Figure 9 in the original paper.

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

%% Extract variables

% Natural log of industrial production
industrial_production_log = data_table.LNIPNSA;

% Natural log of producer prices
producer_price_log = data_table.LNPPINSA;

% Cumulated Romer and Romer monetary policy shock
monetary_shock = data_table.SUMSHCK;

%% Check whether SUMSHCK and SUMSHCKF are identical

same_shock_series = isequal( ...
    data_table.SUMSHCK, ...
    data_table.SUMSHCKF);

fprintf('Are SUMSHCK and SUMSHCKF identical? %d\n', ...
    same_shock_series);

%% Construct the VAR dataset

% Cholesky ordering:
% 1. Industrial production
% 2. Producer price level
% 3. Cumulated monetary policy shock

y1 = [
    industrial_production_log, ...
    producer_price_log, ...
    monetary_shock
    ];

%% Estimate the Cholesky SVAR

lags = 12;

opt1.hor = 48;
opt1.K   = 0;

VAR1 = cvar_(y1, lags, opt1);

%% Extract responses to the monetary policy shock

monetary_shock_index = 3;

% Extract the third structural shock
irfs1 = squeeze( ...
    VAR1.ir_ols(:, :, monetary_shock_index, :) ...
    );

% Ensure that rows represent variables and columns represent horizons
if size(irfs1, 1) ~= 3 && size(irfs1, 2) == 3
    irfs1 = irfs1';
end

% Check the dimensions
if size(irfs1, 1) ~= 3
    error(['The extracted IRF matrix does not have three variables ', ...
        'in its rows. Check the dimensions of VAR1.ir_ols.']);
end

%% Normalise to a one-percentage-point monetary policy innovation

policy_variable_index = 3;

% First column is the impact response at horizon zero
mp_impact_response = irfs1(policy_variable_index, 1);

if abs(mp_impact_response) < 1e-12
    error(['The impact response of SUMSHCK is approximately zero, ', ...
        'so the IRFs cannot be normalised.']);
end

% Desired impact response of the cumulated monetary policy shock
desired_policy_impact = 1;

% Common scaling factor
scaling_factor = ...
    desired_policy_impact / mp_impact_response;

% Apply the same scaling factor to every variable and horizon
monetary_irf_rescaled = ...
    irfs1 * scaling_factor;

%% Convert logged responses into percentages

monetary_irf_for_plot = monetary_irf_rescaled;

% Convert log-point responses into approximate percentage responses
monetary_irf_for_plot(1:2, :) = ...
    100 * monetary_irf_rescaled(1:2, :);

% SUMSHCK remains in percentage points
monetary_irf_for_plot(3, :) = ...
    monetary_irf_rescaled(3, :);

%% Display normalisation information

fprintf('\nOriginal impact response of SUMSHCK: %.6f\n', ...
    mp_impact_response);

fprintf('Scaling factor: %.6f\n', ...
    scaling_factor);

fprintf(['Rescaled impact response of SUMSHCK ', ...
    '(should equal 1): %.6f\n\n'], ...
    monetary_irf_for_plot(3, 1));

%% Construct the horizon vector

number_of_horizons = size(monetary_irf_for_plot, 2);

% The first observation corresponds to horizon zero
horizons = 0:(number_of_horizons - 1);

%% Plot IRFs using the axis limits from the original paper

figure( ...
    'Position', [100, 50, 750, 1050], ...
    'Color', 'w');

plot_layout = tiledlayout(3, 1);

plot_layout.TileSpacing = 'compact';
plot_layout.Padding = 'compact';

%% Panel A: Effect on the cumulated monetary policy shock

ax1 = nexttile;

plot( ...
    ax1, ...
    horizons, ...
    monetary_irf_for_plot(3, :), ...
    'k-', ...
    'LineWidth', 1.5);

hold(ax1, 'on');

yline( ...
    ax1, ...
    0, ...
    'k:', ...
    'LineWidth', 1);

title( ...
    ax1, ...
    'a. Effect on the Cumulated Shock', ...
    'FontSize', 18);

xlabel( ...
    ax1, ...
    'Months after Shock', ...
    'FontSize', 16);

ylabel( ...
    ax1, ...
    'Percentage Points', ...
    'FontSize', 16);

% Original paper's y-axis limits
ylim(ax1, [0, 1.2]);
yticks(ax1, 0:0.2:1.2);

xlim(ax1, [0, 48]);
xticks(ax1, 0:3:48);

ax1.FontSize = 14;
ax1.Box = 'on';
ax1.LineWidth = 1;

%% Panel B: Effect on industrial production

ax2 = nexttile;

plot( ...
    ax2, ...
    horizons, ...
    monetary_irf_for_plot(1, :), ...
    'k-', ...
    'LineWidth', 1.5);

hold(ax2, 'on');

yline( ...
    ax2, ...
    0, ...
    'k:', ...
    'LineWidth', 1);

title( ...
    ax2, ...
    'b. Effect on Output', ...
    'FontSize', 18);

xlabel( ...
    ax2, ...
    'Months after Shock', ...
    'FontSize', 16);

ylabel( ...
    ax2, ...
    'Percent', ...
    'FontSize', 16);

% Original paper's y-axis limits
ylim(ax2, [-5, 2]);
yticks(ax2, -5:1:2);

xlim(ax2, [0, 48]);
xticks(ax2, 0:3:48);

ax2.FontSize = 14;
ax2.Box = 'on';
ax2.LineWidth = 1;

%% Panel C: Effect on the producer price level

ax3 = nexttile;

plot( ...
    ax3, ...
    horizons, ...
    monetary_irf_for_plot(2, :), ...
    'k-', ...
    'LineWidth', 1.5);

hold(ax3, 'on');

yline( ...
    ax3, ...
    0, ...
    'k:', ...
    'LineWidth', 1);

title( ...
    ax3, ...
    'c. Effect on the Price Level', ...
    'FontSize', 18);

xlabel( ...
    ax3, ...
    'Months after Shock', ...
    'FontSize', 16);

ylabel( ...
    ax3, ...
    'Percent', ...
    'FontSize', 16);

% Original paper's y-axis limits
ylim(ax3, [-7, 1]);
yticks(ax3, -7:1:1);

xlim(ax3, [0, 48]);
xticks(ax3, 0:3:48);

ax3.FontSize = 14;
ax3.Box = 'on';
ax3.LineWidth = 1;

%% Ensure all panels use the same x-axis range

linkaxes([ax1, ax2, ax3], 'x');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VAR using monetary shock measure + commodity prices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct the VAR dataset


commodity_prices_log = data_table.LNWCP;

% Make sure the variable name matches the variable created when the
% data were imported. For example:
%
% fed_funds = data_table.FF;

y_MP_COM = [
    industrial_production_log, ...
    producer_price_log, ...
    commodity_prices_log, ...
    monetary_shock
    ];

%% Estimate the Cholesky SVAR

lags = 36;

opt__MP_COM.hor = 48;
opt__MP_COM.K   = 0;

VAR__MP_COM = cvar_(y__MP_COM, lags, opt__MP_COM);

%% Extract responses to the monetary policy shock

% The federal funds rate is the fourth variable, so the fourth
% Cholesky shock is treated as the monetary policy shock.
monetary_shock_index = 4;

irfs__MP_COM = squeeze( ...
    VAR__MP_COM.ir_ols(:, :, monetary_shock_index, :) ...
    );

% Ensure rows represent variables and columns represent horizons
if size(irfs__MP_COM, 1) ~= 4 && size(irfs__MP_COM, 2) == 4
    irfs_ff = irfs__MP_COM';
end

% Check the extracted IRF dimensions
if size(irfs__MP_COM, 1) ~= 4
    error(['The extracted IRF matrix does not have four variables ', ...
        'in its rows. Check the dimensions of VAR_ff.ir_ols.']);
end

%% Normalise to a one-percentage-point federal funds rate innovation

policy_variable_index = 4;

% The first column is the impact response at horizon zero
MP_impact_response = ...
    irfs_MP_COM(policy_variable_index, 1);

if abs(MP_impact_response) < 1e-12
    error(['The impact response of the federal funds rate is ', ...
        'approximately zero, so the IRFs cannot be normalised.']);
end

% Desired impact increase in the federal funds rate
desired_policy_impact = 1;

% Common scaling factor
scaling_factor = ...
    desired_policy_impact / MP_impact_response;

% Apply the same scaling factor to all variables and horizons
irfs_MP_COM_rescaled = ...
    irfs_MP_COM * scaling_factor;

%% Convert logged responses into percentages

irfs_MP_COM_for_plot = irfs_MP_COM_rescaled;

% Industrial production, producer prices and commodity prices are logs.
% Multiply by 100 to express their responses in approximate percent.
irfs_MP_COM_for_plot(1:3, :) = ...
    100 * irfs_MP_COM_rescaled(1:3, :);

% The federal funds rate is already measured in percentage points.
irfs_ff_for_plot(4, :) = ...
    irfs_ff_rescaled(4, :);

%% Display normalisation information

fprintf('\nOriginal impact response of federal funds rate: %.6f\n', ...
    fed_funds_impact_response);

fprintf('Scaling factor: %.6f\n', ...
    scaling_factor);

fprintf(['Rescaled impact response of federal funds rate ', ...
    '(should equal 1): %.6f\n\n'], ...
    irfs_ff_for_plot(4, 1));

%% Construct the horizon vector

number_of_horizons = size(irfs_ff_for_plot, 2);

% First observation corresponds to horizon zero
horizons = 0:(number_of_horizons - 1);

%% Plot the impulse response functions

figure( ...
    'Position', [100, 50, 800, 1250], ...
    'Color', 'w');

plot_layout = tiledlayout(4, 1);

plot_layout.TileSpacing = 'compact';
plot_layout.Padding = 'compact';

%% Panel A: Federal funds rate response

ax1 = nexttile;

plot( ...
    ax1, ...
    horizons, ...
    irfs_ff_for_plot(4, :), ...
    'k-', ...
    'LineWidth', 1.5);

hold(ax1, 'on');

yline( ...
    ax1, ...
    0, ...
    'k:', ...
    'LineWidth', 1);

title( ...
    ax1, ...
    'a. Effect on the Federal Funds Rate', ...
    'FontSize', 18);

xlabel( ...
    ax1, ...
    'Months after Shock', ...
    'FontSize', 16);

ylabel( ...
    ax1, ...
    'Percentage Points', ...
    'FontSize', 16);

ylim(ax1, [0, 1.2]);
yticks(ax1, 0:0.2:1.2);

xlim(ax1, [0, 48]);
xticks(ax1, 0:3:48);

ax1.FontSize = 14;
ax1.Box = 'on';
ax1.LineWidth = 1;

%% Panel B: Industrial production response

ax2 = nexttile;

plot( ...
    ax2, ...
    horizons, ...
    irfs_ff_for_plot(1, :), ...
    'k-', ...
    'LineWidth', 1.5);

hold(ax2, 'on');

yline( ...
    ax2, ...
    0, ...
    'k:', ...
    'LineWidth', 1);

title( ...
    ax2, ...
    'b. Effect on Output', ...
    'FontSize', 18);

xlabel( ...
    ax2, ...
    'Months after Shock', ...
    'FontSize', 16);

ylabel( ...
    ax2, ...
    'Percent', ...
    'FontSize', 16);

ylim(ax2, [-5, 2]);
yticks(ax2, -5:1:2);

xlim(ax2, [0, 48]);
xticks(ax2, 0:3:48);

ax2.FontSize = 14;
ax2.Box = 'on';
ax2.LineWidth = 1;

%% Panel C: Producer price level response

ax3 = nexttile;

plot( ...
    ax3, ...
    horizons, ...
    irfs_ff_for_plot(2, :), ...
    'k-', ...
    'LineWidth', 1.5);

hold(ax3, 'on');

yline( ...
    ax3, ...
    0, ...
    'k:', ...
    'LineWidth', 1);

title( ...
    ax3, ...
    'c. Effect on the Price Level', ...
    'FontSize', 18);

xlabel( ...
    ax3, ...
    'Months after Shock', ...
    'FontSize', 16);

ylabel( ...
    ax3, ...
    'Percent', ...
    'FontSize', 16);

ylim(ax3, [-7, 1]);
yticks(ax3, -7:1:1);

xlim(ax3, [0, 48]);
xticks(ax3, 0:3:48);

ax3.FontSize = 14;
ax3.Box = 'on';
ax3.LineWidth = 1;

%% Panel D: Commodity price response

ax4 = nexttile;

plot( ...
    ax4, ...
    horizons, ...
    irfs_ff_for_plot(3, :), ...
    'k-', ...
    'LineWidth', 1.5);

hold(ax4, 'on');

yline( ...
    ax4, ...
    0, ...
    'k:', ...
    'LineWidth', 1);

title( ...
    ax4, ...
    'd. Effect on Commodity Prices', ...
    'FontSize', 18);

xlabel( ...
    ax4, ...
    'Months after Shock', ...
    'FontSize', 16);

ylabel( ...
    ax4, ...
    'Percent', ...
    'FontSize', 16);

% Commodity-price limits are left automatic because the original
% three-panel figure does not provide limits for this fourth response.
% To set them manually, use, for example:
%
% ylim(ax4, [-10, 5]);
% yticks(ax4, -10:2:5);

xlim(ax4, [0, 48]);
xticks(ax4, 0:3:48);

ax4.FontSize = 14;
ax4.Box = 'on';
ax4.LineWidth = 1;

%% Ensure every panel uses the same x-axis range

linkaxes([ax1, ax2, ax3, ax4], 'x');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Var using fed funds rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Extract variables


Fed_funds = data_table.FF;


%% Construct the VAR dataset

% Cholesky ordering:
% 1. Industrial production
% 2. Producer price level
% 3. Fed Funds

y1 = [
    industrial_production_log, ...
    producer_price_log, ...
    Fed_funds
    ];

%% Estimate the Cholesky SVAR

lags = 36;

opt1.hor = 48;
opt1.K   = 0;

VAR1 = cvar_(y1, lags, opt1);

%% Extract responses to the monetary policy shock

monetary_shock_index = 3;



% Extract the third structural shock
irfs2 = squeeze( ...
    VAR1.ir_ols(:, :, monetary_shock_index, :) ...
    );

% Ensure that rows represent variables and columns represent horizons
if size(irfs2, 1) ~= 3 && size(irfs2, 2) == 3
    irfs2 = irfs2';
end

% Check the dimensions
if size(irfs2, 1) ~= 3
    error(['The extracted IRF matrix does not have three variables ', ...
        'in its rows. Check the dimensions of VAR1.ir_ols.']);
end

%% Normalise to a one-percentage-point monetary policy innovation

policy_variable_index = 3;

% First column is the impact response at horizon zero
mp_impact_response = irfs2(policy_variable_index, 1);

if abs(mp_impact_response) < 1e-12
    error(['The impact response of SUMSHCK is approximately zero, ', ...
        'so the IRFs cannot be normalised.']);
end

% Desired impact response of the cumulated monetary policy shock
desired_policy_impact = 1;

% Common scaling factor
scaling_factor = ...
    desired_policy_impact / mp_impact_response;

% Apply the same scaling factor to every variable and horizon
monetary_irf_rescaled = ...
    irfs2 * scaling_factor;

%% Convert logged responses into percentages

monetary_irf_for_plot = monetary_irf_rescaled;

% Convert log-point responses into approximate percentage responses
monetary_irf_for_plot(1:2, :) = ...
    100 * monetary_irf_rescaled(1:2, :);

% SUMSHCK remains in percentage points
monetary_irf_for_plot(3, :) = ...
    monetary_irf_rescaled(3, :);

%% Display normalisation information

fprintf('\nOriginal impact response of SUMSHCK: %.6f\n', ...
    mp_impact_response);

fprintf('Scaling factor: %.6f\n', ...
    scaling_factor);

fprintf(['Rescaled impact response of SUMSHCK ', ...
    '(should equal 1): %.6f\n\n'], ...
    monetary_irf_for_plot(3, 1));

%% Construct the horizon vector

number_of_horizons = size(monetary_irf_for_plot, 2);

% The first observation corresponds to horizon zero
horizons = 0:(number_of_horizons - 1);

%% Plot IRFs using the axis limits from the original paper

figure( ...
    'Position', [100, 50, 750, 1050], ...
    'Color', 'w');

plot_layout = tiledlayout(3, 1);

plot_layout.TileSpacing = 'compact';
plot_layout.Padding = 'compact';

%% Panel A: Effect on the cumulated monetary policy shock

ax1 = nexttile;

plot( ...
    ax1, ...
    horizons, ...
    monetary_irf_for_plot(3, :), ...
    'k-', ...
    'LineWidth', 1.5);

hold(ax1, 'on');

yline( ...
    ax1, ...
    0, ...
    'k:', ...
    'LineWidth', 1);

title( ...
    ax1, ...
    'a. Effect on the Cumulated Shock', ...
    'FontSize', 18);

xlabel( ...
    ax1, ...
    'Months after Shock', ...
    'FontSize', 16);

ylabel( ...
    ax1, ...
    'Percentage Points', ...
    'FontSize', 16);

% Original paper's y-axis limits
ylim(ax1, [0, 1.2]);
yticks(ax1, 0:0.2:1.2);

xlim(ax1, [0, 48]);
xticks(ax1, 0:3:48);

ax1.FontSize = 14;
ax1.Box = 'on';
ax1.LineWidth = 1;

%% Panel B: Effect on industrial production

ax2 = nexttile;

plot( ...
    ax2, ...
    horizons, ...
    monetary_irf_for_plot(1, :), ...
    'k-', ...
    'LineWidth', 1.5);

hold(ax2, 'on');

yline( ...
    ax2, ...
    0, ...
    'k:', ...
    'LineWidth', 1);

title( ...
    ax2, ...
    'b. Effect on Output', ...
    'FontSize', 18);

xlabel( ...
    ax2, ...
    'Months after Shock', ...
    'FontSize', 16);

ylabel( ...
    ax2, ...
    'Percent', ...
    'FontSize', 16);

% Original paper's y-axis limits
ylim(ax2, [-5, 2]);
yticks(ax2, -5:1:2);

xlim(ax2, [0, 48]);
xticks(ax2, 0:3:48);

ax2.FontSize = 14;
ax2.Box = 'on';
ax2.LineWidth = 1;

%% Panel C: Effect on the producer price level

ax3 = nexttile;

plot( ...
    ax3, ...
    horizons, ...
    monetary_irf_for_plot(2, :), ...
    'k-', ...
    'LineWidth', 1.5);

hold(ax3, 'on');

yline( ...
    ax3, ...
    0, ...
    'k:', ...
    'LineWidth', 1);

title( ...
    ax3, ...
    'c. Effect on the Price Level', ...
    'FontSize', 18);

xlabel( ...
    ax3, ...
    'Months after Shock', ...
    'FontSize', 16);

ylabel( ...
    ax3, ...
    'Percent', ...
    'FontSize', 16);

% Original paper's y-axis limits
ylim(ax3, [-7, 1]);
yticks(ax3, -7:1:1);

xlim(ax3, [0, 48]);
xticks(ax3, 0:3:48);

ax3.FontSize = 14;
ax3.Box = 'on';
ax3.LineWidth = 1;

%% Ensure all panels use the same x-axis range

linkaxes([ax1, ax2, ax3], 'x');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Var using fed funds rate + controlling for commodity prices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Four-variable Cholesky SVAR with the federal funds rate
%
% Cholesky ordering:
% 1. Log industrial production
% 2. Log producer price level
% 3. Log commodity prices
% 4. Federal funds rate
%
% The IRFs are normalised so that the federal funds rate rises by
% one percentage point on impact.

%% Construct the VAR dataset


commodity_prices_log = data_table.LNWCP;

% Make sure the variable name matches the variable created when the
% data were imported. For example:
%
% fed_funds = data_table.FF;

y_ff = [
    industrial_production_log, ...
    producer_price_log, ...
    commodity_prices_log, ...
    Fed_funds
    ];

%% Estimate the Cholesky SVAR

lags = 36;

opt_ff.hor = 48;
opt_ff.K   = 0;

VAR_ff = cvar_(y_ff, lags, opt_ff);

%% Extract responses to the monetary policy shock

% The federal funds rate is the fourth variable, so the fourth
% Cholesky shock is treated as the monetary policy shock.
monetary_shock_index = 4;

irfs_ff = squeeze( ...
    VAR_ff.ir_ols(:, :, monetary_shock_index, :) ...
    );

% Ensure rows represent variables and columns represent horizons
if size(irfs_ff, 1) ~= 4 && size(irfs_ff, 2) == 4
    irfs_ff = irfs_ff';
end

% Check the extracted IRF dimensions
if size(irfs_ff, 1) ~= 4
    error(['The extracted IRF matrix does not have four variables ', ...
        'in its rows. Check the dimensions of VAR_ff.ir_ols.']);
end

%% Normalise to a one-percentage-point federal funds rate innovation

policy_variable_index = 4;

% The first column is the impact response at horizon zero
fed_funds_impact_response = ...
    irfs_ff(policy_variable_index, 1);

if abs(fed_funds_impact_response) < 1e-12
    error(['The impact response of the federal funds rate is ', ...
        'approximately zero, so the IRFs cannot be normalised.']);
end

% Desired impact increase in the federal funds rate
desired_policy_impact = 1;

% Common scaling factor
scaling_factor = ...
    desired_policy_impact / fed_funds_impact_response;

% Apply the same scaling factor to all variables and horizons
irfs_ff_rescaled = ...
    irfs_ff * scaling_factor;

%% Convert logged responses into percentages

irfs_ff_for_plot = irfs_ff_rescaled;

% Industrial production, producer prices and commodity prices are logs.
% Multiply by 100 to express their responses in approximate percent.
irfs_ff_for_plot(1:3, :) = ...
    100 * irfs_ff_rescaled(1:3, :);

% The federal funds rate is already measured in percentage points.
irfs_ff_for_plot(4, :) = ...
    irfs_ff_rescaled(4, :);

%% Display normalisation information

fprintf('\nOriginal impact response of federal funds rate: %.6f\n', ...
    fed_funds_impact_response);

fprintf('Scaling factor: %.6f\n', ...
    scaling_factor);

fprintf(['Rescaled impact response of federal funds rate ', ...
    '(should equal 1): %.6f\n\n'], ...
    irfs_ff_for_plot(4, 1));

%% Construct the horizon vector

number_of_horizons = size(irfs_ff_for_plot, 2);

% First observation corresponds to horizon zero
horizons = 0:(number_of_horizons - 1);

%% Plot the impulse response functions

figure( ...
    'Position', [100, 50, 800, 1250], ...
    'Color', 'w');

plot_layout = tiledlayout(4, 1);

plot_layout.TileSpacing = 'compact';
plot_layout.Padding = 'compact';

%% Panel A: Federal funds rate response

ax1 = nexttile;

plot( ...
    ax1, ...
    horizons, ...
    irfs_ff_for_plot(4, :), ...
    'k-', ...
    'LineWidth', 1.5);

hold(ax1, 'on');

yline( ...
    ax1, ...
    0, ...
    'k:', ...
    'LineWidth', 1);

title( ...
    ax1, ...
    'a. Effect on the Federal Funds Rate', ...
    'FontSize', 18);

xlabel( ...
    ax1, ...
    'Months after Shock', ...
    'FontSize', 16);

ylabel( ...
    ax1, ...
    'Percentage Points', ...
    'FontSize', 16);

ylim(ax1, [0, 1.2]);
yticks(ax1, 0:0.2:1.2);

xlim(ax1, [0, 48]);
xticks(ax1, 0:3:48);

ax1.FontSize = 14;
ax1.Box = 'on';
ax1.LineWidth = 1;

%% Panel B: Industrial production response

ax2 = nexttile;

plot( ...
    ax2, ...
    horizons, ...
    irfs_ff_for_plot(1, :), ...
    'k-', ...
    'LineWidth', 1.5);

hold(ax2, 'on');

yline( ...
    ax2, ...
    0, ...
    'k:', ...
    'LineWidth', 1);

title( ...
    ax2, ...
    'b. Effect on Output', ...
    'FontSize', 18);

xlabel( ...
    ax2, ...
    'Months after Shock', ...
    'FontSize', 16);

ylabel( ...
    ax2, ...
    'Percent', ...
    'FontSize', 16);

ylim(ax2, [-5, 2]);
yticks(ax2, -5:1:2);

xlim(ax2, [0, 48]);
xticks(ax2, 0:3:48);

ax2.FontSize = 14;
ax2.Box = 'on';
ax2.LineWidth = 1;

%% Panel C: Producer price level response

ax3 = nexttile;

plot( ...
    ax3, ...
    horizons, ...
    irfs_ff_for_plot(2, :), ...
    'k-', ...
    'LineWidth', 1.5);

hold(ax3, 'on');

yline( ...
    ax3, ...
    0, ...
    'k:', ...
    'LineWidth', 1);

title( ...
    ax3, ...
    'c. Effect on the Price Level', ...
    'FontSize', 18);

xlabel( ...
    ax3, ...
    'Months after Shock', ...
    'FontSize', 16);

ylabel( ...
    ax3, ...
    'Percent', ...
    'FontSize', 16);

ylim(ax3, [-7, 1]);
yticks(ax3, -7:1:1);

xlim(ax3, [0, 48]);
xticks(ax3, 0:3:48);

ax3.FontSize = 14;
ax3.Box = 'on';
ax3.LineWidth = 1;

%% Panel D: Commodity price response

ax4 = nexttile;

plot( ...
    ax4, ...
    horizons, ...
    irfs_ff_for_plot(3, :), ...
    'k-', ...
    'LineWidth', 1.5);

hold(ax4, 'on');

yline( ...
    ax4, ...
    0, ...
    'k:', ...
    'LineWidth', 1);

title( ...
    ax4, ...
    'd. Effect on Commodity Prices', ...
    'FontSize', 18);

xlabel( ...
    ax4, ...
    'Months after Shock', ...
    'FontSize', 16);

ylabel( ...
    ax4, ...
    'Percent', ...
    'FontSize', 16);

% Commodity-price limits are left automatic because the original
% three-panel figure does not provide limits for this fourth response.
% To set them manually, use, for example:
%
% ylim(ax4, [-10, 5]);
% yticks(ax4, -10:2:5);

xlim(ax4, [0, 48]);
xticks(ax4, 0:3:48);

ax4.FontSize = 14;
ax4.Box = 'on';
ax4.LineWidth = 1;

%% Ensure every panel uses the same x-axis range

linkaxes([ax1, ax2, ax3, ax4], 'x');