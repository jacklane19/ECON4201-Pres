%% Cholesky SVAR using the Romer and Romer monetary policy shock
% The IRFs are normalised to a one-percentage-point innovation in SUMSHCK.
% Responses of logged output and prices are converted into percentages.
%
% The VAR includes:
%   - A constant, included automatically by cvar_
%   - Eleven monthly seasonal dummies
%   - January as the omitted reference month
%
% IMPORTANT:
% The exogenous-variable section inside cvar_.m must align the exogenous
% variables using:
%
% exogenous(idx(1) + lags : idx(end), :)
%
% rather than:
%
% exogenous(idx, :)

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

number_of_observations = height(data_table);

%% Construct the monthly date series

variable_names = data_table.Properties.VariableNames;

if ismember('DATE', variable_names)

    raw_dates = data_table.DATE;

    if isdatetime(raw_dates)

        monthly_dates = raw_dates;

    elseif isnumeric(raw_dates)

        finite_dates = raw_dates(isfinite(raw_dates));

        % Check whether DATE is stored as YYYYMM, such as 196901
        appears_to_be_yyyymm = ...
            ~isempty(finite_dates) && ...
            all(finite_dates >= 180001) && ...
            all(finite_dates <= 220012) && ...
            all(mod(finite_dates, 1) == 0);

        if appears_to_be_yyyymm

            year_number = floor(raw_dates / 100);
            month_number_from_date = mod(raw_dates, 100);

            monthly_dates = datetime( ...
                year_number, ...
                month_number_from_date, ...
                1);

        else

            % Otherwise, treat the values as Excel serial dates
            monthly_dates = datetime( ...
                raw_dates, ...
                'ConvertFrom', 'excel');

        end

    else

        % Convert character, cell or string dates
        monthly_dates = datetime(string(raw_dates));

    end

elseif all(ismember({'YEAR', 'MONTH'}, variable_names))

    monthly_dates = datetime( ...
        data_table.YEAR, ...
        data_table.MONTH, ...
        1);

else

    error([ ...
        'The data table does not contain a recognised DATE variable ', ...
        'or separate YEAR and MONTH variables.']);

end

monthly_dates = monthly_dates(:);

if length(monthly_dates) ~= number_of_observations

    error(['The date series contains %d observations, but the data ', ...
        'table contains %d observations.'], ...
        length(monthly_dates), ...
        number_of_observations);

end

if any(isnat(monthly_dates))

    error('The monthly date series contains invalid or missing dates.');

end

%% Extract the endogenous variables

% Natural log of industrial production
industrial_production_log = data_table.LNIPNSA(:);

% Natural log of the producer price index
producer_price_log = data_table.LNPPINSA(:);

% Cumulated Romer and Romer monetary policy shock
monetary_shock = data_table.SUMSHCK(:);

% Original Romer and Romer monetary policy shock
original_rr_shock = data_table.RESID(:);

% Additional variables, not used in this specification
commodity_prices_log = data_table.LNWCP(:);
fed_funds = data_table.FF(:);

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
%
% Placing the monetary policy shock last means that output and prices
% cannot respond contemporaneously to the monetary policy shock.

y1 = [
    industrial_production_log, ...
    producer_price_log, ...
    monetary_shock
    ];

%% Check the endogenous data

if any(~isfinite(y1), 'all')

    error('The VAR dataset contains missing or non-finite observations.');

end

%% Construct eleven monthly seasonal dummies

% cvar_ includes a constant automatically.
%
% January is omitted to avoid perfect multicollinearity:
%
% Column 1  = February
% Column 2  = March
% ...
% Column 11 = December

calendar_month = month(monthly_dates);

seasonal_dummies = zeros(number_of_observations, 11);

for current_month = 2:12

    seasonal_dummies(:, current_month - 1) = ...
        double(calendar_month == current_month);

end

seasonal_dummy_names = {
    'February'
    'March'
    'April'
    'May'
    'June'
    'July'
    'August'
    'September'
    'October'
    'November'
    'December'
    };

fprintf(['Included %d monthly seasonal dummies ', ...
    '(February to December).\n'], ...
    size(seasonal_dummies, 2));

fprintf('January is the omitted reference month.\n');

%% Check the seasonal dummy matrix

if size(seasonal_dummies, 1) ~= size(y1, 1)

    error(['The seasonal dummy matrix contains %d rows, while y1 ', ...
        'contains %d rows.'], ...
        size(seasonal_dummies, 1), ...
        size(y1, 1));

end

if size(seasonal_dummies, 2) ~= 11

    error('The seasonal dummy matrix should contain 11 columns.');

end

if any(~isfinite(seasonal_dummies), 'all')

    error('The seasonal dummy matrix contains invalid observations.');

end

%% Estimate the Cholesky SVAR

lags = 36;

opt1.hor = 48;

% No bootstrap simulations
opt1.K = 0;

% Add the eleven seasonal dummies to every VAR equation.
% Do not add these dummies directly to y1.
opt1.exogenous = seasonal_dummies;

fprintf('\nEstimating VAR with:\n');
fprintf('  %d endogenous variables\n', size(y1, 2));
fprintf('  %d monthly lags\n', lags);
fprintf('  %d seasonal dummies\n', size(seasonal_dummies, 2));
fprintf('  %d total observations\n\n', size(y1, 1));

VAR1 = cvar_(y1, lags, opt1);

%% Extract the IRFs to the monetary policy shock

monetary_shock_index = 3;

% After squeeze, the expected dimensions are:
% variable x horizon
irfs1 = squeeze( ...
    VAR1.ir_ols(:, :, monetary_shock_index, :) ...
    );

% Ensure variables are stored in rows
if size(irfs1, 1) ~= size(y1, 2) && ...
        size(irfs1, 2) == size(y1, 2)

    irfs1 = irfs1';

end

if size(irfs1, 1) ~= size(y1, 2)

    error(['The extracted IRF matrix has unexpected dimensions: ', ...
        '%d rows by %d columns.'], ...
        size(irfs1, 1), ...
        size(irfs1, 2));

end

%% Normalise to a one-percentage-point innovation in SUMSHCK

policy_variable_index = 3;

% The first IRF column is the impact response at horizon zero
mp_impact_response = irfs1(policy_variable_index, 1);

if abs(mp_impact_response) < 1e-12

    error(['The impact response of the monetary policy variable is ', ...
        'approximately zero, so the IRFs cannot be normalised.']);

end

% Desired impact response of SUMSHCK
desired_policy_impact = 1;

% Common normalisation factor
scaling_factor = ...
    desired_policy_impact / mp_impact_response;

% Apply the same scaling factor to every variable and horizon
monetary_irf_rescaled = ...
    irfs1 * scaling_factor;

%% Convert logged-variable responses into percentages

monetary_irf_for_plot = monetary_irf_rescaled;

% LNIPNSA and LNPPINSA are natural logarithms.
% Multiplying by 100 converts log-point responses into approximate
% percentage responses.
monetary_irf_for_plot(1:2, :) = ...
    100 * monetary_irf_rescaled(1:2, :);

% SUMSHCK is already measured in percentage points
monetary_irf_for_plot(3, :) = ...
    monetary_irf_rescaled(3, :);

%% Verify the IRF normalisation

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

%% Extract the residual from the monetary-policy equation

% VAR1.e_ols contains the reduced-form residuals from each VAR equation.
% Column 3 is the residual from the SUMSHCK equation.

monetary_shock_VAR_residuals = VAR1.e_ols(:, 3);
monetary_shock_VAR_residuals = ...
    monetary_shock_VAR_residuals(:);

%% Align the original Romer-Romer shock with the VAR residuals

% The first 36 observations are lost because the VAR includes 36 lags.

monetary_shock_RR = original_rr_shock(lags + 1:end);
monetary_shock_RR = monetary_shock_RR(:);

residual_dates = monthly_dates(lags + 1:end);
residual_dates = residual_dates(:);

%% Check residual-series dimensions

fprintf('Romer-Romer shock observations: %d\n', ...
    length(monetary_shock_RR));

fprintf('VAR residual observations: %d\n', ...
    length(monetary_shock_VAR_residuals));

fprintf('Residual date observations: %d\n', ...
    length(residual_dates));

if length(monetary_shock_RR) ~= ...
        length(monetary_shock_VAR_residuals)

    error(['The shock series have different lengths: ', ...
        'Romer-Romer shock = %d, VAR residual = %d.'], ...
        length(monetary_shock_RR), ...
        length(monetary_shock_VAR_residuals));

end

if length(residual_dates) ~= ...
        length(monetary_shock_VAR_residuals)

    error(['The residual date series has %d observations, while ', ...
        'the VAR residual series has %d observations.'], ...
        length(residual_dates), ...
        length(monetary_shock_VAR_residuals));

end

%% Remove missing observations for the correlation calculation

valid_residual_observations = ...
    isfinite(monetary_shock_RR) & ...
    isfinite(monetary_shock_VAR_residuals);

if sum(valid_residual_observations) < 3

    error(['There are not enough valid observations to calculate ', ...
        'the correlation.']);

end

rr_shock_valid = ...
    monetary_shock_RR(valid_residual_observations);

var_residual_valid = ...
    monetary_shock_VAR_residuals(valid_residual_observations);

dates_valid = ...
    residual_dates(valid_residual_observations);

%% Calculate the residual correlation

[residual_correlation, correlation_p_value] = corr( ...
    rr_shock_valid, ...
    var_residual_valid);

fprintf('\nCorrelation between the original Romer-Romer shock\n');
fprintf('and the VAR monetary-policy residual: %.4f\n', ...
    residual_correlation);

fprintf('Correlation p-value: %.6f\n\n', ...
    correlation_p_value);

%% Plot the two residual series over time

figure;

plot( ...
    dates_valid, ...
    rr_shock_valid, ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Original Romer-Romer shock');

hold on;

plot( ...
    dates_valid, ...
    var_residual_valid, ...
    'LineWidth', 1.2, ...
    'DisplayName', 'VAR monetary-policy residual');

yline( ...
    0, ...
    '--', ...
    'LineWidth', 1, ...
    'HandleVisibility', 'off');

hold off;

title( ...
    sprintf(['Original Romer-Romer Shock and VAR Residual ', ...
    '(r = %.4f)'], residual_correlation), ...
    'FontSize', 18);

xlabel( ...
    'Date', ...
    'FontSize', 16);

ylabel( ...
    'Monetary policy shock', ...
    'FontSize', 16);

legend( ...
    'Location', 'best', ...
    'FontSize', 12);

set(gca, ...
    'FontSize', 14);

grid on;
box on;

%% Plot the correlation as a scatter plot

% The horizontal axis contains the original Romer-Romer residual.
% The vertical axis contains the residual from the third VAR equation.

figure;

scatter( ...
    rr_shock_valid, ...
    var_residual_valid, ...
    30, ...
    'filled', ...
    'DisplayName', 'Monthly observations');

hold on;

% Add an ordinary least-squares fitted line
fitted_line_coefficients = polyfit( ...
    rr_shock_valid, ...
    var_residual_valid, ...
    1);

fitted_line_x = linspace( ...
    min(rr_shock_valid), ...
    max(rr_shock_valid), ...
    100);

fitted_line_y = polyval( ...
    fitted_line_coefficients, ...
    fitted_line_x);

plot( ...
    fitted_line_x, ...
    fitted_line_y, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Linear fitted line');

hold off;

title( ...
    sprintf(['Correlation Between Original Shock and VAR Residual ', ...
    '(r = %.4f, p = %.4f)'], ...
    residual_correlation, ...
    correlation_p_value), ...
    'FontSize', 18);

xlabel( ...
    'Original Romer-Romer shock (RESID)', ...
    'FontSize', 16);

ylabel( ...
    'Residual from the VAR SUMSHCK equation', ...
    'FontSize', 16);

legend( ...
    'Location', 'best', ...
    'FontSize', 12);

set(gca, ...
    'FontSize', 14);

grid on;
box on;