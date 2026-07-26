clear;
clc;


% Setting up path for empirical macro toolbox
addpath(fullfile(pwd, '..', 'BVAR_', 'bvartools'));



%% Read CSV file as a MATLAB table

data_table = readtable("RomerandRomerDataAppendix.xls",'Sheet',"DATA BY MONTH");

% Display the first few rows
head(data_table)

% Display the imported variable names
disp(data_table.Properties.VariableNames)

%
industrial_production_log = data_table.LNIPNSA;
producer_price_log = data_table.LNPPINSA;
commodity_prices_log = data_table.LNWCP;
monetary_shock = data_table.SUMSHCK;
fed_funds = data_table.FF;


 isequal(data_table.SUMSHCK, data_table.SUMSHCKF)
 %Series are not the same, but based off the descriptions I'm note sure
 %what the difference is





%% Cholesky Decomposition


%Settings applying to each 
%-36 lags
%-48 month Horizon for IRF
%-No confidence intervals (but may add at some point just to see how
%sensative

lags = 36;


%First specification

opt1.hor = 48;
opt1.K = 0;

y1 = [industrial_production_log, producer_price_log, monetary_shock];

VAR1 = cvar_(y1, lags, opt1);

irfs1 = squeeze(VAR1.ir_ols(:, :, 3, :));

opt1_plot.varnames     = {'industrial_production_log','producer_price_log', 'monetary_shock'};
opt1_plot.shocksnames  = {'monetary_shock'};
opt1_plot.nplots = [1, 2];

plot_irfs_(irfs1, opt1_plot);




%%
%Second specification




%Third Specification


lags = 4;

opt1.hor = 20;
opt1.K = 1;
opt2.hor = 20;
opt2.K=1;

y1 = [output, inflation, interestRate];
y2 = [output, consumption, investment, wages, hours, inflation, interestRate];



VAR1 = cvar_(y1, lags, opt1);
VAR2 = cvar_(y2, lags, opt2);

irfs1 = VAR1.ir_ols([1, 2], :, 3, :);
irfs2 = VAR2.ir_ols([1, 6], :, 7, :);

% Plot IRF for first system
opt1_plot.varnames     = {'Output', 'Inflation'};
opt1_plot.shocksnames  = {'Interest Rate'};
opt1_plot.nplots = [1, 2];

plot_irfs_(irfs1, opt1_plot);

% Adjust axes for comparability
fig = gcf;
ax1  = findall(fig, 'Type', 'axes');
ax1  = flipud(ax1(:));      % findall returns axes in reverse creation order

ylim(ax1(1), [output_y_axis(2), output_y_axis(1)]);      % Output subplot
ylim(ax1(2), [inf_y_axis(2), inf_y_axis(1)]);   % Inflation subplot




% Plot IRF for second system
opt2_plot = opt1_plot;
opt2_plot.varnames     = {'Output', 'Inflation'};

plot_irfs_(irfs2, opt2_plot);

% Adjust axes for comparability
fig = gcf;
ax1  = findall(fig, 'Type', 'axes');
ax1  = flipud(ax1(:));      % findall returns axes in reverse creation order

ylim(ax1(1), [output_y_axis(2), output_y_axis(1)]);      % Output subplot
ylim(ax1(2), [inf_y_axis(2), inf_y_axis(1)]);   % Inflation subplot

% Save
saveas(fig, './Q3_plots/Q3 IRFs System 2.png');

% Get Frobenius norm
Sigma_u = VAR1.Sigma_ols;
P = chol(Sigma_u, 'lower');
frob = norm(Sigma_u- P*P', 'fro');
disp(frob);



%Tests of validity

