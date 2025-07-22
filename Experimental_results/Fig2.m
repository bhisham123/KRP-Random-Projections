clc
clear all
close all

load("Cauchy_relative_error_250_2.mat")

figure;
plot(ranks(1:end), mean(err_rsthosvd_krp(1:end,:),2), '--r', 'LineWidth', 2, 'DisplayName', 'RSTHOSVD-KRP'); hold on
plot(ranks(1:end), mean(err_rhosvd_krp_memo(1:end,:),2), '-r','LineWidth', 2, 'DisplayName', 'RHOSVD-KRP-MEMO'); 
plot(ranks(1:end), mean(err_rhosvd_krp(1:end,:),2), '-','Color', '#ff7f0e', 'LineWidth', 2, 'DisplayName', 'RHOSVD-KRP'); 
plot(ranks(1:end), mean(err_rsthosvd(1:end,:),2),'--b', 'LineWidth', 2, 'DisplayName', 'RSTHOSVD');
plot(ranks(1:end), mean(err_rhosvd(1:end,:),2), '-b', 'LineWidth', 2, 'DisplayName', 'RHOSVD');
plot(ranks(1:end), mean(err_sthosvd(1:end,:),2), '--k', 'LineWidth', 2, 'DisplayName', 'STHOSVD');
plot(ranks(1:end), mean(err_hosvd(1:end,:),2), '-k','LineWidth', 2, 'DisplayName', 'HOSVD');


% Add labels, legend, title
legend('Location', 'northeast');
xlabel('Rank');
xticks([0,50,100,150])
ylabel('Relative Error');
set(gca, 'YScale', 'log'); 
set(gca, 'FontSize', 16);
title('Accuracy')
ylim([1e-16 max(err_rhosvd_krp_memo, [],'all')*1.5]);   % Set y-axis from 0 to 10

print(gcf, 'Fig_2.png', '-dpng', '-r300')
% print(gcf, 'Fig_2.pdf', '-dpng', '-r300')
