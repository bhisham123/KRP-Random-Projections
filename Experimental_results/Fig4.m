clc
clear all
close all

load tangaroa_pqr.mat 


f = figure;
f.Position =  [100, 100, 1450, 500]; 
subplot(1, 2, 1)
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_hosvd_gram)./cellfun(@(x) sum(mean(x(2:end,:), 1)),time_rhosvd_krp_memo), '-or', 'LineWidth', 2.5, 'DisplayName', 'RHOSVD-KRP-MEMO'); hold on
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_hosvd_gram)./cellfun(@(x) sum(mean(x(2:end,:), 1)),time_rhosvd_krp), '-*','Color', '#FFA500', 'LineWidth', 2.5, 'DisplayName', 'RHOSVD-KRP'); 
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_hosvd_gram)./cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rhosvd), '-x','Color','#0074D9', 'LineWidth', 2.5, 'DisplayName', 'RHOSVD');
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_hosvd_gram)./cellfun(@(x) sum(mean(x(2:end,:), 1)), time_svdsketch), '-+', 'Color', '#8c564b','LineWidth', 2.5, 'DisplayName', 'RHOSVD-svdsketch');
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_hosvd_gram)./cellfun(@(x) sum(mean(x(2:end,:), 1)), time_hosvd_gram), '-sk', 'LineWidth', 2.5, 'DisplayName', 'HOSVD');
% Adding labels and legend
legend('Location', 'northeast','FontSize',15);
xlabel('Rank')
xticks(ranks)
ylabel('Speedup (w.r.t. HOSVD)');
set(gca, 'FontSize', 16.5);
title('Compression of training data')

subplot(1, 2, 2)
plot(ranks, cellfun(@(x) mean(x(:)),err_rhosvd_krp_mttkrps), '-or', 'LineWidth', 2.5, 'DisplayName', 'RHOSVD-KRP-MEMO'); hold on;
plot(ranks, cellfun(@(x) mean(x(:)),err_rhosvd_krp), '-*','Color', '#FFA500','LineWidth', 2.5, 'DisplayName', 'RHOSVD-KRP'); 
plot(ranks, cellfun(@(x) mean(x(:)), err_rhosvd), '-x', 'Color', '#0074D9','LineWidth', 2.5, 'DisplayName', 'RHOSVD');
plot(ranks, cellfun(@(x) mean(x(:)), err_svdsketch), '-+','Color', '#8c564b', 'LineWidth', 2.5, 'DisplayName', 'RHOSVD-svdsketch');
plot(ranks, cellfun(@(x) mean(x(:)),err_hosvd_gram), '-sk', 'LineWidth', 2.5, 'DisplayName', 'HOSVD');

% Add labels, legend, title
legend('Location', 'northeast','FontSize',15);
xlabel('Rank');
xticks(ranks)
ylabel('Relative Error');
% set(gca, 'YScale', 'log'); 
title('Test error');
set(gca, 'FontSize', 16.5);

file_name = "Fig4.png";
print(gcf, file_name, '-dpng', '-r300')
