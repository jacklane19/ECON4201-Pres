%% Cholesky SVAR using the Romer and Romer monetary policy shock
% The IRFs are normalised to a one-percentage-point innovation in SUMSHCK.
% Responses of logged output and prices are converted into percentages.

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

%% Extract the variables

% Natural log of industrial production
industrial_production_log = data_table.LNIPNSA;

% Natural log of the producer price index
producer_price_log = data_table.LNPPINSA;

% Cumulated Romer and Romer monetary policy shock
monetary_shock = data_table.SUMSHCK;

% Additional variables, not used in this specification
commodity_prices_log = data_table.LNWCP;
fed_funds = data_table.FF;

% Check whether SUMSHCK and SUMSHCKF are identical
same_shock_series = isequal(data_table.SUMSHCK, data_table.SUMSHCKF);

fprintf('Are SUMSHCK and SUMSHCKF identical? %d\n', ...
    same_shock_series);

%% Construct the VAR dataset

% Cholesky ordering:
% 1. Industrial production
% 2. Producer price level
% 3. Monetary policy shock
%
% Placing the monetary policy shock last means that output and prices
% cannot respond contemporaneously to the monetary policy shock.

y1 = [
    industrial_production_log, ...
    producer_price_log, ...
    monetary_shock
    ];

%% Estimate the Cholesky SVAR

lags = 36;

opt1.hor = 48;
opt1.K   = 0;

VAR1 = cvar_(y1, lags, opt1);

%% Extract the IRFs to the monetary policy shock

monetary_shock_index = 3;

% After squeeze, the expected dimensions are:
% variable x horizon
irfs1 = squeeze(VAR1.ir_ols(:, :, monetary_shock_index, :));

% Ensure variables are stored in rows
if size(irfs1, 1) ~= size(y1, 2) && ...
        size(irfs1, 2) == size(y1, 2)

    irfs1 = irfs1';

end

%% Normalise to a one-percentage-point innovation in SUMSHCK

policy_variable_index = 3;

% The first IRF column is the impact response at horizon zero
mp_impact_response = irfs1(policy_variable_index, 1);

% Check that the impact response is nonzero
if abs(mp_impact_response) < 1e-12
    error(['The impact response of the monetary policy variable is ', ...
        'approximately zero, so the IRFs cannot be normalised.']);
end

% Desired impact response of SUMSHCK
desired_policy_impact = 1;

% Common normalisation factor
scaling_factor = desired_policy_impact / mp_impact_response;

% Apply the same scaling factor to every variable and horizon
monetary_irf_rescaled = irfs1 * scaling_factor;

%% Convert logged-variable responses into percentages

monetary_irf_for_plot = monetary_irf_rescaled;

% LNIPNSA and LNPPINSA are natural logarithms.
% Multiplying their IRFs by 100 converts log-point responses into
% approximate percentage responses.
monetary_irf_for_plot(1:2, :) = ...
    100 * monetary_irf_rescaled(1:2, :);

% Do not multiply SUMSHCK by 100 because it is already measured in
% percentage points.
monetary_irf_for_plot(3, :) = ...
    monetary_irf_rescaled(3, :);

%% Verify the normalisation

fprintf('\nOriginal impact response of SUMSHCK: %.6f\n', ...
    mp_impact_response);

fprintf('Scaling factor: %.6f\n', ...
    scaling_factor);

fprintf(['Rescaled impact response of SUMSHCK ', ...
    '(should equal 1): %.6f\n\n'], ...
    monetary_irf_for_plot(3, 1));

%% Plot the rescaled IRFs

opt1_plot.varnames = {
    'Industrial production'
    'Producer price level'
    'Cumulated monetary shock'
    };

opt1_plot.shocksnames = {
    'Monetary policy shock'
    };

opt1_plot.nplots = [1, 3];

plot_irfs_(monetary_irf_for_plot, opt1_plot);