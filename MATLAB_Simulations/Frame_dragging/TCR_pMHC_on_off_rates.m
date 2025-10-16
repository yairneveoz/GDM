% TCR_pMHC_on_off_rates

clc
clear 

log10_on_rates = linspace(-2,0,21);
log10_off_rates = linspace(-6,0,20);

[LOG10_on_rates,LOG10_off_rates] = ...
    meshgrid(log10_on_rates,log10_off_rates);

% contact_duration = 10.^LOG10_on_rates.*10.^LOG10_off_rates;
contact_duration = LOG10_on_rates - LOG10_off_rates;

figure(2)
surf(log10_on_rates,log10_off_rates,contact_duration,...
    'EdgeColor','none')
colormap(hot)
cb = colorbar();
% cb.Ticks = [-2:2:8];
cb.Ticks = [-2:2:6];
cb.TickLabels = ({10.^cb.Ticks});
title(cb,'Sec','FontSize',12,'Rotation',0)
axis square
view(2)
xlabel('P_{on} (sec^{-1})')
ylabel('P_{off} (sec^{-1})')
title('Duration time of TCR-pMHC contact','FontSize',12)
Xticks = [-2, -1, 0];
xticks(Xticks)
% yticks([-9, -6, -3, 0])
xticklabels({10.^Xticks})
% yticklabels({'10^{-9}', '10^{-6}', '10^{-3}', '1'})
Yticks = [-6, -4, -2, 0];
yticks(Yticks)
% yticklabels({'10^{-6}', '10^{-4}', '10^{-2}', '1'})
yticklabels({10.^Yticks})

