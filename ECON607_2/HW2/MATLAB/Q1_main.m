%% ECON607_II - HW2 - Q1
%
%  Anirudh Yadav
%  18 March 2018
%
%  NOTES: The hours measure is constructed from BLS Productivity and Costs survey.
%  I take the raw index measure (retrieved from FRED) and deflate by the
%  non-institutional civilian population. I tried to use Shimer's preferred
%  measure (from CPS data), but I couldn't find a seasonally adjusted
%  version.
%
%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

clear;
close all;


%----------------------------------------------------------------
% 1. Import and prepare data
%----------------------------------------------------------------

% Import raw data
names = {"gdp", "c", "i", "hours_1","w", "tfp"};
for ind = 1:length(names)
  data.(names{ind}) = xlsread(strcat('/Users/Anirudh/Desktop/PhD/ECON607/HW2/Data/', names{ind}, '.xlsx'));
end

% Log the data
for ind = 1:length(names)
  logdata.(names{ind}) = log(data.(names{ind})(:,2));
end

% Create a time variable
t = datetime(1948,1,1):calmonths(3):datetime(2017,10,1);

%----------------------------------------------------------------
% 2. Compute trend data and detrended data using HPFILTER
%----------------------------------------------------------------

for ind = 1:length(names)
  trend_data.(names{ind}) = hpfilter(logdata.(names{ind}),1600);
  detrended_data.(names{ind}) = logdata.(names{ind})-trend_data.(names{ind});
end

%----------------------------------------------------------------
% 3. Plot data
%----------------------------------------------------------------

% Plot raw log data and HP trend
for ind = 1:length(names)
    plot_1 = figure
    plot(t,logdata.(names{ind}),t,trend_data.(names{ind}));
    title(strcat("log ",names{ind}),'Interpreter', 'none')
    legend('raw','HP trend','Location','northwest')
    saveas(plot_1, strcat("plot1_",names{ind},".eps"))
end

% Plot detrended data
for ind = 1:length(names)
    plot_2 = figure
    plot(t,detrended_data.(names{ind}));
    title(strcat("Detrended ",names{ind}),'Interpreter', 'none')
    saveas(plot_2, strcat("plot2_",names{ind},".eps"))
end

%----------------------------------------------------------------
% 4. Compute standard deviations of detrended data
%----------------------------------------------------------------

% Compute stdevs for each series
for ind = 1:length(names)
  stdev.(names{ind}) = std(detrended_data.(names{ind}));
end

% Compute stdevs relative to gdp
for ind = 1:length(names)
  stdevrel.(names{ind}) = stdev.(names{ind})/stdev.gdp;
end

%----------------------------------------------------------------
% 5. Compute first-order autocorrelations
%----------------------------------------------------------------

% Create lagged detrended data
for ind = 1:length(names)
  lagged_data_1.(names{ind}) = lagmatrix(detrended_data.(names{ind}),1);
end

% Compute first-order autocorrelations
for ind = 1:length(names)
    corr_lag1.(names{ind})=corrcoef(detrended_data.(names{ind}),lagged_data_1.(names{ind}),'rows','pairwise');
end

%----------------------------------------------------------------
% 6. Compute contemporaneous correlations with gdp
%----------------------------------------------------------------

for ind = 1:length(names)
    corr0_Y.(names{ind})=corrcoef(detrended_data.gdp,detrended_data.(names{ind}),'rows','pairwise');
end