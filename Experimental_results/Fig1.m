clc
clear all
close all

load("SysID_results.mat")


% ---- Compute average Hausdorff distances ----
Time_rsvd = Time_rsvd(:,2:end);
Time_gauss = Time_gauss(:,2:end);
Time_krp = Time_krp(:,2:end);

avg_rsvd = mean(Time_rsvd,2);
avg_gauss = mean(Time_gauss,2);
avg_krp = mean(Time_krp,2);


% ---- Compute average Hausdorff distances ----
hd_rsvd = hd_rsvd(:,2:end);
hd_gauss = hd_gauss(:,2:end);
hd_krp = hd_krp(:,2:end);

avg_hd_rsvd  = mean(hd_rsvd, 2);
avg_hd_gauss = mean(hd_gauss, 2);
avg_hd_krp   = mean(hd_krp, 2);



figure;
set(gcf, 'Position', [100, 100, 1000, 350]);  % [left, bottom, width, height]
% ---- Subplot 1: Hausdorff Distance ----
subplot(1,2,2)
plot(Ranks, abs(avg_hd_krp),   'r*-', 'DisplayName', 'R-KRP', 'LineWidth', 2); hold on 
plot(Ranks, abs(avg_hd_gauss), 'ks-', 'DisplayName', 'R-Gauss', 'LineWidth', 2); 
plot(Ranks, abs(avg_hd_rsvd),   'bo-', 'DisplayName', 'RandERA', 'LineWidth', 2); 

set(gca, 'FontSize', 14.8);
xlabel('Rank (r)');
ylabel('Average Hausdorff Distance');
legend('show', 'Location', 'northeast');
title('Error')
xticks(Ranks);


% ---- Subplot 2: Average Runtime ----
subplot(1,2,1)
plot(Ranks, avg_krp,    'r*-', 'DisplayName', 'R-KRP', 'LineWidth', 2); hold on
plot(Ranks, avg_gauss,  'ks-', 'DisplayName', 'R-Gauss', 'LineWidth', 2); 
plot(Ranks, avg_rsvd,   'bo-', 'DisplayName', 'RandERA', 'LineWidth', 2); 
set(gca, 'FontSize', 14.8);
xlabel('Rank (r)');
ylabel('Average Runtime (Seconds)');
legend('show', 'Location', 'northwest');
xticks(Ranks);
title('Runtime')
file_name = "Fig1.png";
% print(gcf, file_name, '-dpng', '-r300');

disp("Speeedup over R-Gauss")
avg_gauss./avg_krp
disp("Speedup over R-RandERA")
avg_rsvd./avg_krp
