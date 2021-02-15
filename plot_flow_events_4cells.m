clear all
close all
clc

load('flow_events_4cells.mat')

kcoh = [0 0.003 0.03 0.3];
num_runs = 25;

for i = 1:4
    flow_rev_mean(i) = mean(flow_reversals(i,:));
    flow_rev_std(i) = std(flow_reversals(i,:));
    flow_drop_mean(i) = mean(flow_drop(i,:));
    flow_drop_std(i) = std(flow_drop(i,:));
end

figure(1), hold on, box on, grid on
errorbar(log10(kcoh(2:4)), flow_rev_mean(2:4), 1.96*flow_rev_std(2:4)/sqrt(num_runs), 'b.-', 'MarkerSize', 20, 'LineWidth', 3)
errorbar(log10(kcoh(2:4)), flow_drop_mean(2:4), 1.96*flow_drop_std(2:4)/sqrt(num_runs), 'r.-', 'MarkerSize', 20, 'LineWidth', 3)

errorbar(-3, flow_rev_mean(1), 1.96*flow_rev_std(1)/sqrt(num_runs), 'bx', 'MarkerSize', 10, 'LineWidth', 3)
errorbar(-3, flow_drop_mean(1), 1.96*flow_drop_std(1)/sqrt(num_runs), 'rx', 'MarkerSize', 10, 'LineWidth', 3)

axis([-3 0 0 1000])
set(gca, 'FontSize', 24)
set(gca, 'LineWidth', 2)
xlabel(' log_{10}(k_{coh}) ', 'FontSize', 24)
ylabel(' events ', 'FontSize', 24)
set(gcf, 'Color', 'w')

ax = gca;
ax.XTick = [-2.5 -1.5 -0.5];
ax.YTick = [0 250 500 750 1000];

legend('flow reversals (w/ cohesion)', 'flow loss (w/ cohesion)', 'flow reversals (w/o cohesion)', 'flow loss (w/o cohesion)', 'location', 'northeastoutside')
legend boxoff

fig = gcf;
pos = fig.Position;
fig.Position = [1 2 1280 640];

flow_rev_no_cohesion = flow_reversals(1,:);
flow_drop_no_cohesion = flow_drop(1,:);

disp('flow reversals')
[h, p] = ttest2(flow_rev_no_cohesion, flow_reversals(2,:), 'VarType', 'unequal')
[h, p] = ttest2(flow_rev_no_cohesion, flow_reversals(3,:), 'VarType', 'unequal')
[h, p] = ttest2(flow_rev_no_cohesion, flow_reversals(4,:), 'VarType', 'unequal')

disp('flow drops')
[h, p] = ttest2(flow_drop_no_cohesion, flow_drop(2,:), 'VarType', 'unequal')
[h, p] = ttest2(flow_drop_no_cohesion, flow_drop(3,:), 'VarType', 'unequal')
[h, p] = ttest2(flow_drop_no_cohesion, flow_drop(4,:), 'VarType', 'unequal')

