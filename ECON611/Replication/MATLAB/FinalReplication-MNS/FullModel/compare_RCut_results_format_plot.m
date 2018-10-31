title(titletext,'FontSize',Format.FontSize)
L = legend('Incomplete Markets', 'Complete Markets', 'Location','Best');
set(L,'FontSize',Format.FontSize);
xlabel('Quarter','FontSize',Format.FontSize);
if Format.interestrate
    ylabel('Basis Points','FontSize',Format.FontSize);
else
    ylabel('Percentage Points','FontSize',Format.FontSize);
end

for i = 1:2
    set(p(i),'color',Format.colors{i});
    set(p(i),'linestyle',Format.styles{i});
    set(p(i),'linewidth',Format.widths{i});
end