clc
clear all
close all

load("Cauchy_time_250_2.mat")

% Create figure and set size
f = figure;
f.Position = [100, 100, 1460, 550]; 
subplot(1, 2, 1)
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rsthosvd_krp), '--dr', 'LineWidth', 2, 'DisplayName', 'RSTHOSVD-KRP'); hold on;
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)),time_rhosvd_krp_memo), '-+r', 'LineWidth', 2, 'DisplayName', 'RHOSVD-KRP-MEMO'); 
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)),time_rhosvd_krp), '-*','Color', '#ff7f0e', 'LineWidth', 2, 'DisplayName', 'RHOSVD-KRP'); 
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rsthosvd), '--pb', 'LineWidth', 2, 'DisplayName', 'RSTHOSVD');
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rhosvd), '-ob', 'LineWidth', 2, 'DisplayName', 'RHOSVD');
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_sthosvd_new),  '--xk', 'LineWidth', 2, 'DisplayName', 'STHOSVD');
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_hosvd_new), '-sk', 'LineWidth', 2, 'DisplayName', 'HOSVD');

% Adding labels and legend
xlabel('Rank');
xticks(ranks)
ylabel('Average Runtime (Seconds)');
ylim([min(cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rsthosvd_krp))-10, max(cellfun(@(x) sum(mean(x(2:end,:), 1)), time_hosvd))+50])
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 16);
lgd = legend('Location', 'southeast');
lgd.FontSize = 15.1;
% grid on;
title('Runtime')


subplot(1, 2, 2)
num_ranks = numel(ranks);
data_all = [];
x_labels = {};
methods = {'RSTHOSVD-KRP','RHOSVD-KRP-MEMO','RHOSVD-KRP', 'RSTHOSVD', 'RHOSVD'};

for i = 1:num_ranks
    % Average the time data for each method (excluding first element if it's for initialization)
    t_rhosvd_krp = mean(time_rhosvd_krp{i}(2:end,:), 1);
    t_rhosvd_krp_memo = mean(time_rhosvd_krp_memo{i}(2:end,:), 1);
    t_rhosvd = mean(time_rhosvd{i}(2:end,:), 1);
    t_rsthosvd_krp = mean(time_rsthosvd_krp{i}(2:end,:), 1);
    t_rsthosvd = mean(time_rsthosvd{i}(2:end,:), 1);

    % Append data for each method
    data_all = [data_all; t_rsthosvd_krp; t_rhosvd_krp_memo; t_rhosvd_krp; t_rsthosvd; t_rhosvd];

    % Append method labels for the x-axis
    for j = 1:numel(methods)
        x_labels{end+1} = sprintf('%s (r=%d)', methods{j}, ranks(i));
    end
     % Insert a blank (zero) row to create visual space
    if i < num_ranks
        data_all = [data_all; zeros(1, size(t_rsthosvd_krp, 2))];
        x_labels{end+1} = '';  % Blank label for spacing
    end
end

data_all_new = data_all;
data_all_new(:,4) = data_all(:,5);
data_all_new(:,5) = data_all(:,4);
% Plot the stacked bar chart
bar(data_all_new, 'stacked');
% Label the y-axis
ylabel('Time (Seconds)');
% Set x-axis labels and rotate them for better readability
xticks(1:numel(x_labels));  % Set positions of x-ticks to match the number of labels
xticklabels(x_labels);       % Set the labels at the x-tick positions
xtickangle(45);              % Rotate labels for readability
set(gca, 'FontSize', 16); 
lgd = legend('Core','mult/MTTKRP', 'QR', 'Unfolding','RNG',  'Location', 'northwest');
lgd.FontSize = 15;
ax = gca;
ax.XAxis.FontSize = 12.5;
title('Time Breakdown')

file_name = "Fig3.png";
saveas(gcf, file_name);


% %% individual plots
% clear all
% close all
% 
% load("Cauchy_time_250_2.mat")
% 
% % Create figure and set size
% f = figure;
% plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rsthosvd_krp), '--dr', 'LineWidth', 2, 'DisplayName', 'RSTHOSVD-KRP'); hold on;
% plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)),time_rhosvd_krp_memo), '-+r', 'LineWidth', 2, 'DisplayName', 'RHOSVD-KRP-MEMO'); 
% plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)),time_rhosvd_krp), '-*','Color', '#ff7f0e', 'LineWidth', 2, 'DisplayName', 'RHOSVD-KRP'); 
% plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rsthosvd), '--pb', 'LineWidth', 2, 'DisplayName', 'RSTHOSVD');
% plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rhosvd), '-ob', 'LineWidth', 2, 'DisplayName', 'RHOSVD');
% plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_sthosvd_new),  '--xk', 'LineWidth', 2, 'DisplayName', 'STHOSVD');
% plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_hosvd_new), '-sk', 'LineWidth', 2, 'DisplayName', 'HOSVD');
% 
% % Adding labels and legend
% xlabel('Rank');
% xticks(ranks)
% ylabel('Average Runtime (Seconds)');
% ylim([min(cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rsthosvd_krp))-10, max(cellfun(@(x) sum(mean(x(2:end,:), 1)), time_hosvd))+50])
% set(gca, 'YScale', 'log');
% set(gca, 'FontSize', 14);
% lgd = legend('Location', 'southeast','FontSize',14);
% % lgd.FontSize = 15;
% title('Runtime')
% % print(gcf, 'Runtime_Cauchy_2.png', '-dpng', '-r300')
% 
% 
% 
% figure('Position', [560  420  560  500]);  % 800x600 px figure at (100, 100)
% num_ranks = numel(ranks);
% % Initialize data storage and labels
% data_all = [];
% x_labels = {};
% methods = {'RSTHOSVD-KRP','RHOSVD-KRP-MEMO','RHOSVD-KRP', 'RSTHOSVD', 'RHOSVD'};
% 
% % Loop through ranks and methods
% for i = 1:num_ranks
%     % Average the time data for each method (excluding first element if it's for initialization)
%     t_rhosvd_krp = mean(time_rhosvd_krp{i}(2:end,:), 1);
%     t_rhosvd_krp_memo = mean(time_rhosvd_krp_memo{i}(2:end,:), 1);
%     t_rhosvd = mean(time_rhosvd{i}(2:end,:), 1);
%     t_rsthosvd_krp = mean(time_rsthosvd_krp{i}(2:end,:), 1);
%     t_rsthosvd = mean(time_rsthosvd{i}(2:end,:), 1);
% 
%     % Append data for each method
%     data_all = [data_all; t_rsthosvd_krp; t_rhosvd_krp_memo; t_rhosvd_krp; t_rsthosvd; t_rhosvd];
% 
%     % Append method labels for the x-axis
%     for j = 1:numel(methods)
%         x_labels{end+1} = sprintf('%s (r=%d)', methods{j}, ranks(i));
%     end
%      % Insert a blank (zero) row to create visual space
%     if i < num_ranks
%         data_all = [data_all; zeros(1, size(t_rsthosvd_krp, 2))];
%         x_labels{end+1} = '';  % Blank label for spacing
%     end
% end
% 
% data_all_new = data_all;
% data_all_new(:,4) = data_all(:,5);
% data_all_new(:,5) = data_all(:,4);
% % Plot the stacked bar chart
% bar(data_all_new, 'stacked');
% % Label the y-axis
% ylabel('Time (Seconds)');
% % Set x-axis labels and rotate them for better readability
% xticks(1:numel(x_labels));  % Set positions of x-ticks to match the number of labels
% xticklabels(x_labels);       % Set the labels at the x-tick positions
% xtickangle(55);              % Rotate labels for readability
% set(gca, 'FontSize', 14); 
% lgd = legend('Core','mult/MTTKRP', 'QR', 'Unfolding','RNG',  'Location', 'northwest');
% % lgd.FontSize = 14;
% ax = gca;
% ax.XAxis.FontSize = 11.5;
% title('Time Breakdown')
% % 
% % file_name = "Running_Time_split_p="+num2str(l)+"_250.png";
% % print(gcf, 'Time_Breakdown_Cauchy_2.png', '-dpng', '-r300')
% 
% 



