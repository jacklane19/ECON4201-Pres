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







monetary_shock_RR =  data_table.RESID(lags + 1:end);


Monetary_shock_VAR_residuals = VAR1.e_ols(:, 3);

% Remove the first 36 observations from the original series
% because they are lost when estimating the VAR

correlation = corr(monetary_shock_RR, Monetary_shock_VAR_residuals);

fprintf('Correlation: %.4f\n', correlation);

size(monetary_shock_RR)
size(Monetary_shock_VAR_residuals)

%% Extract and align the two monetary-policy shock series

% Romer and Romer residual series
% Remove the first observations lost due to the VAR lags
monetary_shock_RR = data_table.RESID(lags + 1:end);

% Residuals from the monetary-policy equation of the VAR
Monetary_shock_VAR_residuals = VAR1.e_ols(:, 3);

% Ensure both series are column vectors
monetary_shock_RR = monetary_shock_RR(:);
Monetary_shock_VAR_residuals = Monetary_shock_VAR_residuals(:);

%% Check that the series have equal lengths

fprintf('Romer-Romer shock observations: %d\n', ...
    length(monetary_shock_RR));

fprintf('VAR residual observations: %d\n', ...
    length(Monetary_shock_VAR_residuals));

if length(monetary_shock_RR) ~= length(Monetary_shock_VAR_residuals)
    error(['The series have different lengths: RR = %d, ' ...
        'VAR residuals = %d.'], ...
        length(monetary_shock_RR), ...
        length(Monetary_shock_VAR_residuals));
end

%% Calculate the correlation

correlation = corr( ...
    monetary_shock_RR, ...
    Monetary_shock_VAR_residuals, ...
    'Rows', 'complete');

fprintf('Correlation: %.4f\n', correlation);

%% Plot both series on the same chart

% Map residuals to their observations in the original dataset
observation_index = ...
    (lags + 1):(lags + length(monetary_shock_RR));

figure;

plot(observation_index, monetary_shock_RR, ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Romer and Romer shock');

hold on;

plot(observation_index, Monetary_shock_VAR_residuals, ...
    'LineWidth', 1.2, ...
    'DisplayName', 'VAR monetary-policy residual');

yline(0, '--', ...
    'LineWidth', 1, ...
    'DisplayName', 'Zero line');

hold off;

title('Comparison of Monetary Policy Shock Series', ...
    'FontSize', 18);

xlabel('Observation in Original Dataset', ...
    'FontSize', 16);

ylabel('Monetary Policy Shock', ...
    'FontSize', 16);

legend('Location', 'best');

set(gca, 'FontSize', 14);

grid on;
box on;


