%% Quandt-Andrews test for the three-variable Romer-Romer VAR


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

% Cholesky ordering:
% 1. Output level
% 2. Price level
% 3. Monetary-policy measure

% Replace these variable names if your workspace uses different names.
y = [
    industrial_production_log, ...
    producer_price_log, ...
    monetary_shock
    ];

% Romer and Romer monthly VAR lag length
lags = 36;

% Use actual monthly dates where available.
%
% Example:
% dates = datetime(data_table.YEAR, data_table.MONTH, 1);
%
% When no date variable has been constructed, observation numbers can be
% used temporarily:
dates = (1:size(y, 1))';

results_QA = quandtAndrewsVAR( ...
    y, ...
    lags, ...
    dates, ...
    'Trim', 0.15, ...
    'MinExtraObservations', 5, ...
    'NumSimulations', 4999, ...
    'Seed', 4201, ...
    'Plot', true);

%% Display main results

disp(results_QA.Summary);
disp(results_QA.CriticalValues);
