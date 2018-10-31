% compare R-cut transition results
close all; clear all; clc;


%---FOMATTING-----

PRESENTATION = true;

if PRESENTATION
    Format.colors = {'red', 'black' };
    Format.styles = {'-','-'};
    Format.widths = {2,2};
    Format.figsize = [460 200 700 525];
    Format.FontSize = 16;
else
    gray = repmat(0.6,1,3);
    Format.colors = {[0 0 0], [0 0 0]};
    Format.styles = {'-','--'};
    Format.widths = {1,1};
    Format.figsize = [460 200 560 420];
    Format.FontSize = 12;
end




load('Results/IMOutputGE.mat');
load('Results/CMOutputGE.mat');

load('Results/IMInflationGE.mat');
load('Results/CMInflationGE.mat');

plot_T = 40;

%---- Y ------
fighandle = figure('Position',Format.figsize);
axes('FontSize',Format.FontSize);
p = plot(0:plot_T,100*[ IMOutputGE(1:plot_T+1)' CMOutputGE(1:plot_T+1)']);
ylim([-0.05 0.3])

titletext = 'Output';
    
Format.interestrate = false;
compare_RCut_results_format_plot;
    
    
if PRESENTATION
    saveas(fighandle,['Figures/RCut' titletext  '_Presentation.pdf']);
else
    saveas(fighandle,['Figures/RCut' titletext  '_Paper.pdf']);
end



%---- inflation ------
fighandle = figure('Position',Format.figsize);
axes('FontSize',Format.FontSize);
p = plot(0:plot_T,100*[ IMInflationGE(1:plot_T+1)' CMInflationGE(1:plot_T+1)']);

titletext = 'Inflation';
    
Format.interestrate = false;
compare_RCut_results_format_plot;
    
    
if PRESENTATION
    saveas(fighandle,['Figures/RCut' titletext  '_Presentation.pdf']);
else
    saveas(fighandle,['Figures/RCut' titletext  '_Paper.pdf']);
end





