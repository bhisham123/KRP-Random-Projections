clc 
clear all
maxNumCompThreads(1);

addpath("./tensor_toolbox/")
addpath("./tucker/")

% % Generate a 3D array with size 4x4x4 using l = 2 (Euclidean norm)
alpha = 2
mode_size = 250;
n_modes = 4;
T = tensor(generate_pnorm_array(n_modes, mode_size, alpha));

%% set params
% ranks = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150]; %ranks for accuracy plot
ranks = [10,20,50,100,150]; %ranks for time plot
p = 0; %oversampling
runs = 5; %indepdent runs
modes=1:4;

%% initialize timings
time_rsthosvd = cell(length(ranks),1);
time_rsthosvd_krp = cell(length(ranks),1);
time_sthosvd = cell(length(ranks),1);


time_rhosvd = cell(length(ranks),1);
time_rhosvd_krp = cell(length(ranks),1);
time_rhosvd_krp_memo = cell(length(ranks),1);
time_hosvd = cell(length(ranks),1);


err_rsthosvd = zeros(length(ranks),runs);
err_rsthosvd_krp = zeros(length(ranks),runs);
err_sthosvd = zeros(length(ranks),runs);


err_rhosvd = zeros(length(ranks),runs);
err_rhosvd_krp = zeros(length(ranks),runs);
err_rhosvd_krp_memo = zeros(length(ranks),runs);
err_hosvd = zeros(length(ranks),runs);



normT = norm(T);

for j = 1:length(ranks)
    r = [ranks(j),ranks(j),ranks(j),ranks(j)];

    time_rsthosvd{j} = zeros(runs,5);
    time_rsthosvd_krp{j} = zeros(runs,5);
    time_sthosvd{j} = zeros(runs,5);
    time_rhosvd{j} = zeros(runs,5);
    time_rhosvd_krp{j} = zeros(runs,5);
    time_rhosvd_krp_memo{j} = zeros(runs,5);
    time_hosvd{j} = zeros(runs,5);
    
    disp(['RSTHOSVD-KRP (Rank = ', num2str(ranks(j)),')'])
    for i = 1:runs
        [T1,times] = rsthosvd_krp(T,r,p,modes);
        time_rsthosvd_krp{j}(i,:) = times;
        err_rsthosvd_krp(j,i) =  norm(full(T1)-T)/normT;
    end
    
    disp(['RSTHOSVD (Rank = ', num2str(ranks(j)),')'])
    for i = 1:runs
        [T1,times] = rsthosvd(T,r,p,modes);
        time_rsthosvd{j}(i,:) = times;
        err_rsthosvd(j,i) =  norm(full(T1)-T)/normT;
    end

    disp(['STHOSVD (Rank = ', num2str(ranks(j)),')'])
    for i = 1:runs
        [T1,times] = sthosvd(T,r,modes);
        time_sthosvd{j}(i,:) = times;
        err_sthosvd(j,i) = norm(full(T1)-T)/normT;
    end


    disp(['RHOSVD-KRP (Rank = ', num2str(ranks(j)),')'])
    for i = 1:runs
        [T1,times] = rhosvd_krp(T,r,p,modes);
        time_rhosvd_krp{j}(i,:) = times;
        err_rhosvd_krp(j,i) =  norm(full(T1)-T)/normT;
    end

    disp(['RHOSVD-KRP-MEMOIZATION (Rank = ', num2str(ranks(j)),')'])
    for i = 1:runs
        [T1,times] = rhosvd_krp_mttkrps(T,r,p,modes);
        time_rhosvd_krp_memo{j}(i,:) = times;
        err_rhosvd_krp_memo(j,i) =  norm(full(T1)-T)/normT;
    end
    
    disp(['RHOSVD (Rank = ', num2str(ranks(j)),')'])
    for i = 1:runs
        [T1,times] = rhosvd(T,r,p,modes);
        time_rhosvd{j}(i,:) = times;
        err_rhosvd(j,i) =  norm(full(T1)-T)/normT;
    end

    disp(['HOSVD (Rank = ', num2str(ranks(j)),')'])
    for i = 1:runs
        [T1,times] = hosvd1(T,r,modes);
        time_hosvd{j}(i,:) = times;
        err_hosvd(j,i) = norm(full(T1)-T)/normT;
    end
end

file_name = "./Experimental_results/Cauchy_"+num2str(mode_size)+"_"+ num2str(alpha)+".mat"; 
save(file_name,"time_rhosvd_krp_memo","time_rhosvd_krp","time_rhosvd","time_hosvd","time_rsthosvd_krp","time_rsthosvd","time_sthosvd","err_hosvd","err_rhosvd","err_rhosvd_krp","err_rhosvd_krp_memo","err_rsthosvd_krp","err_rsthosvd","err_sthosvd","ranks","alpha",'mode_size')

%% Create figure and set size
f = figure;
f.Position = [100, 100, 1460, 550]; 
subplot(1, 2, 1)
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rsthosvd_krp), '--dr', 'LineWidth', 2, 'DisplayName', 'RSTHOSVD-KRP'); hold on;
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)),time_rhosvd_krp_memo), '-+r', 'LineWidth', 2, 'DisplayName', 'RHOSVD-KRP-MEMO'); 
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)),time_rhosvd_krp), '-*','Color', '#ff7f0e', 'LineWidth', 2, 'DisplayName', 'RHOSVD-KRP'); 
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rsthosvd), '--pb', 'LineWidth', 2, 'DisplayName', 'RSTHOSVD');
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rhosvd), '-ob', 'LineWidth', 2, 'DisplayName', 'RHOSVD');
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_sthosvd),  '--xk', 'LineWidth', 2, 'DisplayName', 'STHOSVD');
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_hosvd), '-sk', 'LineWidth', 2, 'DisplayName', 'HOSVD');

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
% Initialize data storage and labels
data_all = [];
x_labels = {};
methods = {'RSTHOSVD-KRP','RHOSVD-KRP-MEMO','RHOSVD-KRP', 'RSTHOSVD', 'RHOSVD'};

% Loop through ranks and methods
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
% 
% file_name = "Running_Time_split_"+num2str(mode_size)+"_"+num2str(alpha)+".png";
% print(gcf, file_name,'-dpng', '-r300')





% Individual plots of running time and time split
f = figure;
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rsthosvd_krp), '--dr', 'LineWidth', 2, 'DisplayName', 'RSTHOSVD-KRP'); hold on;
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)),time_rhosvd_krp_memo), '-+r', 'LineWidth', 2, 'DisplayName', 'RHOSVD-KRP-MEMO'); 
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)),time_rhosvd_krp), '-*','Color', '#ff7f0e', 'LineWidth', 2, 'DisplayName', 'RHOSVD-KRP'); 
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rsthosvd), '--pb', 'LineWidth', 2, 'DisplayName', 'RSTHOSVD');
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rhosvd), '-ob', 'LineWidth', 2, 'DisplayName', 'RHOSVD');
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_sthosvd),  '--xk', 'LineWidth', 2, 'DisplayName', 'STHOSVD');
plot(ranks, cellfun(@(x) sum(mean(x(2:end,:), 1)), time_hosvd), '-sk', 'LineWidth', 2, 'DisplayName', 'HOSVD');
% Adding labels and legend
xlabel('Rank');
xticks(ranks)
ylabel('Average Runtime (Seconds)');
ylim([min(cellfun(@(x) sum(mean(x(2:end,:), 1)), time_rsthosvd_krp))-10, max(cellfun(@(x) sum(mean(x(2:end,:), 1)), time_hosvd))+50])
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);
lgd = legend('Location', 'southeast');
% lgd.FontSize = 15;
title('Runtime')
file_name = "Runtime_Cauchy_"+num2str(mode_size)+"_"+num2str(alpha)+".png";
print(gcf, file_name, '-dpng', '-r300')



figure;
num_ranks = numel(ranks);
% Initialize data storage and labels
data_all = [];
x_labels = {};
methods = {'RSTHOSVD-KRP','RHOSVD-KRP-MEMO','RHOSVD-KRP', 'RSTHOSVD', 'RHOSVD'};

% Loop through ranks and methods
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
xtickangle(55);              % Rotate labels for readability
set(gca, 'FontSize', 14); 
lgd = legend('Core','mult/MTTKRP', 'QR', 'Unfolding','RNG',  'Location', 'northwest');
% lgd.FontSize = 14;
ax = gca;
ax.XAxis.FontSize = 11.5;
title('Time Breakdown')
file_name = "Time_Breakdown_Cauchy_"+num2str(mode_size)+"_"+num2str(alpha)+".png";
print(gcf, file_name, '-dpng', '-r300')



% Relative error plot
f = figure;
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
set(gca, 'FontSize', 14);
title('Accuracy')

% file_name = "Relative_Error_Cauchy_"+num2str(mode_size)+"_"+num2str(alpha)+".png";
% print(gcf, file_name, '-dpng', '-r300')





function X = generate_pnorm_array(n, size_per_dim, p)
    % Validate p
    if p < 1
        error('p must be greater than or equal to 1');
    end

    % Generate a cell array of index grids
    indices = cell(1, n);
    [indices{:}] = ndgrid(1:size_per_dim);

    % Initialize the output array
    X = zeros(size_per_dim * ones(1, n));

    % Compute the p-norm at each index
    norm_p = zeros(size(X));
    for k = 1:n
        norm_p = norm_p + abs(indices{k}).^p;
    end
    norm_p = norm_p.^(1/p);

    % Avoid division by zero (though indices start from 1, so norm_p should never be zero)
    X = 1 ./ norm_p;
end
