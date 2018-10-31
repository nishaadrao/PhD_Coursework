% compare ZLB transition results
close all; clear all; clc;


%---FOMATTING-----

PRESENTATION = true;

if PRESENTATION
    Format.colors = {'red', 'red','black', 'black', };
    Format.styles = {'-','--','-','--'};
    Format.widths = {2,2,2,2};
    Format.figsize = [460 200 700 525];
    Format.FontSize = 16;
else
    gray = repmat(0.6,1,3);
    Format.colors = {[0 0 0], [0 0 0], gray, gray};
    Format.styles = {'-','--','-','--'};
    Format.widths = {1,1,1,1};
    Format.figsize = [460 200 560 420];
    Format.FontSize = 12;
end




IncmpMkts.naive = load('Results/transition_ZLB_naive.mat');
IncmpMkts.extended = load('Results/transition_ZLB_extended.mat');
CmpMkts.naive = load('Results/transition_ZLB_cmp_mkts_naive.mat');
CmpMkts.extended = load('Results/transition_ZLB_cmp_mkts_extended.mat');

plot_T = 40;
relative = @(x)(100*x(2:plot_T+1)/x(1)-100);


%---- Y ------
fighandle = figure('Position',Format.figsize);
axes('FontSize',Format.FontSize);
p = plot([  relative(IncmpMkts.naive.Y)' ...
        relative(IncmpMkts.extended.Y)' ...
        relative(CmpMkts.naive.Y)' ...
        relative(CmpMkts.extended.Y)']);

titletext = 'Output';
    
Format.interestrate = false;
compare_ZLB_results_format_plot;
    
    
if PRESENTATION
    saveas(fighandle,['Figures/ZLB' titletext  '_Presentation.pdf']);
else
    saveas(fighandle,['Figures/ZLB' titletext  '_Paper.pdf']);
end

tmp = relative(IncmpMkts.naive.Y)' - relative(CmpMkts.naive.Y)';
tmp(1)

%---- inflation ------
fighandle = figure('Position',Format.figsize);
axes('FontSize',Format.FontSize);
p = plot([  relative(IncmpMkts.naive.ppi)' ...
        relative(IncmpMkts.extended.ppi)' ...
        relative(CmpMkts.naive.ppi)' ...
        relative(CmpMkts.extended.ppi)']);

titletext = 'Inflation';
    
Format.interestrate = false;
compare_ZLB_results_format_plot;
    
    
if PRESENTATION
    saveas(fighandle,['Figures/ZLB' titletext  '_Presentation.pdf']);
else
    saveas(fighandle,['Figures/ZLB' titletext  '_Paper.pdf']);
end



%---- nominal interest rate ------
fighandle = figure('Position',Format.figsize);
axes('FontSize',Format.FontSize);
p = plot(10000*[  IncmpMkts.naive.ii(2:plot_T)' ...
        IncmpMkts.extended.ii(2:plot_T)' ...
        CmpMkts.naive.ii(2:plot_T)' ...
        CmpMkts.extended.ii(2:plot_T)']);

xlim([15 35])
    
titletext = 'Nominal Interest Rate';
    
Format.interestrate = true;
compare_ZLB_results_format_plot;

    
titletext(ismember(titletext,' ,.:;!')) = [];    
if PRESENTATION
    saveas(fighandle,['Figures/ZLB' titletext  '_Presentation.pdf']);
else
    saveas(fighandle,['Figures/ZLB' titletext  '_Paper.pdf']);
end




