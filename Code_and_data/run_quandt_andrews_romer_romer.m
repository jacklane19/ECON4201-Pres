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
lags = 12;

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

%% Plot the individual VAR series for visual break inspection

% Variables:
% 1. Output level
% 2. Price level
% 3. Monetary-policy measure

series_data = {
    industrial_production_log
    producer_price_log
    monetary_shock
    };

series_titles = {
    'Output Level'
    'Price Level'
    'Monetary Policy Measure'
    };

y_axis_labels = {
    'Log industrial production'
    'Log producer price level'
    'Romer-Romer monetary policy measure'
    };

%% Construct the date vector

% Use your actual monthly dates when available.
%
% Example:
% dates = datetime(data_table.YEAR, data_table.MONTH, 1);

% Use observation numbers if no date vector has been constructed:
dates = (1:length(industrial_production_log))';

%% Plot settings

moving_average_window = 12;
show_moving_average   = true;
show_estimated_break  = true;

% Determine whether an estimated breakpoint is available
break_available = ...
    exist('results_QA', 'var') && ...
    isfield(results_QA, 'EstimatedBreakDate') && ...
    ~isempty(results_QA.EstimatedBreakDate);

%% Create one figure for each series

for series_number = 1:length(series_data)

    current_series = series_data{series_number};

    % Ensure both variables are column vectors
    current_series = current_series(:);
    plot_dates      = dates(:);

    % Remove observations missing from this particular series
    valid_observations = ...
        ~isnan(current_series) & ...
        ~ismissing(plot_dates);

    current_series = current_series(valid_observations);
    current_dates  = plot_dates(valid_observations);

    %% Create figure

    figure( ...
        'Position', ...
        [100, 100, 1250, 650]);

    % Plot the original monthly series
    plot( ...
        current_dates, ...
        current_series, ...
        'LineWidth', ...
        1.2, ...
        'DisplayName', ...
        'Monthly series');

    hold on;

    %% Add a 12-month moving average

    if show_moving_average

        moving_average = movmean( ...
            current_series, ...
            moving_average_window, ...
            'omitnan');

        plot( ...
            current_dates, ...
            moving_average, ...
            'LineWidth', ...
            2, ...
            'DisplayName', ...
            '12-month moving average');

    end

    %% Mark the estimated Quandt-Andrews breakpoint

    if show_estimated_break && break_available

        estimated_break_date = ...
            results_QA.EstimatedBreakDate;

        xline( ...
            estimated_break_date, ...
            '--', ...
            'LineWidth', ...
            1.5, ...
            'DisplayName', ...
            'Estimated VAR breakpoint');

    end

    %% Figure formatting

    title( ...
        [series_titles{series_number}, ...
         ': Visual Inspection for Structural Breaks'], ...
        'FontSize', ...
        18);

    xlabel( ...
        'Date', ...
        'FontSize', ...
        16);

    ylabel( ...
        y_axis_labels{series_number}, ...
        'FontSize', ...
        16);

    current_axes = gca;
    current_axes.FontSize = 14;

    xlim([current_dates(1), current_dates(end)]);

    grid on;
    box on;

    legend( ...
        'Location', ...
        'best', ...
        'FontSize', ...
        14);

    hold off;

end