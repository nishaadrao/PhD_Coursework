close all; clear all; clc;


% FORMATTING ---

PRESENTATION = false;


if PRESENTATION
    Format.FontSize = 16;
    Format.colors = {'k','r'};
    Format.widths = {2,2};
    Format.styles = {'-','-'};
    Format.figsize = [460 200 700 525];
else
    Format.colors = {[0 0 0], [0 0 0]};
    Format.styles = {'--','-'};
    Format.widths = {1,1};
    Format.figsize = [460 200 560 420];
    Format.FontSize = 12;
end



load Results/HorizonTable_IM_Working.mat
HT_IM = HorizonTable;

load Results/HorizonTable_CM_Working.mat
HT_CM = HorizonTable;

I = HT_CM(:,1) <= 40;
HT_CM = HT_CM(I,:);
HT_CM(:,2:3) =  100*(HT_CM(:,2:3)-1);

I = HT_CM(:,3) > 2;
HT_CM(I,3) = NaN;

I = HT_IM(:,1) <= 40 & HT_IM(:,1) > 0;
HT_IM = HT_IM(I,:);
HT_IM(:,2:3) =  100*(HT_IM(:,2:3)-1);

fh = figure('Position',Format.figsize);
axes1 = axes('Parent',fh,'FontSize',Format.FontSize);
p = zeros(1,2);
p(1) = plot(HT_CM(:,1),HT_CM(:,2));
hold on
p(2) = plot(HT_IM(:,1),HT_IM(:,2));
hold off
ylim([0 .30])
xlim([0 40])

for i = 1:2
    set(p(i),'LineStyle',Format.styles{i})
    set(p(i),'color',Format.colors{i})
    set(p(i),'linewidth',Format.widths{i})
end
ylabel('Percentage Points','FontSize',Format.FontSize)
xlabel('Horizon in Quarters','FontSize',Format.FontSize)
title('Output','FontSize',Format.FontSize)
L = legend('Complete Markets', 'Incomplete Markets', 'location','Best');
set(L,'FontSize',Format.FontSize);
    
if PRESENTATION
    saveas(fh,['Figures/MultiHorizonOutput_Presentation.pdf']);
else
    saveas(fh,['Figures/MultiHorizonOutput_Paper.pdf']);
end




fh = figure('Position',Format.figsize);
axes2 = axes('Parent',fh,'FontSize',14);
p = zeros(1,2);
p(1) = plot(HT_CM(:,1),HT_CM(:,3));
hold on
p(2) = plot(HT_IM(:,1),HT_IM(:,3));
hold off
xlim([0 40])


for i = 1:2
    set(p(i),'LineStyle',Format.styles{i})
    set(p(i),'color',Format.colors{i})
    set(p(i),'linewidth',Format.widths{i})
end
ylabel('Percentage Points','FontSize',Format.FontSize)
xlabel('Horizon in Quarters','FontSize',Format.FontSize)
title('Inflation','FontSize',Format.FontSize)
L = legend('Complete Markets','Incomplete Markets',  'location','Best');
set(L,'FontSize',Format.FontSize);

if PRESENTATION
    saveas(fh,['Figures/MultiHorizonInflation_Presentation.pdf']);
else
    saveas(fh,['Figures/MultiHorizonInflation_Paper.pdf']);
end