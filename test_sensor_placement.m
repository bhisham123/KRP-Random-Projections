% simulation data from https://cgl.ethz.ch/research/visualization/data.php 
clear
maxNumCompThreads(1);
addpath('./tensor_toolbox/')
addpath("./sensor_placement/")
%% sim velocity data on 2D grid over time
udata = ncread('tangaroa.nc','u');  %Data from https://cgl.ethz.ch/Downloads/Data/ScientificData/tangaroa3d_nc.zip
% size(udata) % 300 x 180 x 120 x 201

%% Split into training and testing data
uctrain = double(udata(1:2:end, 1:2:end, 1:2:end,1:150));
uctest = double(udata(1:2:end, 1:2:end, 1:2:end,151:end));
%% Perform mean subtraction
ucm = mean(uctrain, 4);
uctrain = bsxfun(@minus, uctrain, ucm);
clear udata
%% convert to tensor objects 
uctrain = tensor(uctrain);
uctest = tensor(uctest);
size(uctrain) %150    90    60   150
%% set params
ranks = [5,10,15,20,25,30]; %ranks
p = 0; %oversampling
runs = 10; %indepdent runs
modes=1:3;
indSelectionMethod = 'pqr'; %subset selection method
%% initialize timings 
time_rhosvd_krp_memo  = cell(1,length(ranks));
time_rhosvd_krp = cell(1,length(ranks));
time_rhosvd = cell(1,length(ranks));
time_hosvd_gram  = cell(1,length(ranks));
time_svdsketch  = cell(1,length(ranks));

%% To store test errror
err_rhosvd_krp_mttkrps  = cell(length(ranks),1);
err_rhosvd_krp  = cell(length(ranks),1);
err_rhosvd = cell(length(ranks),1);
err_hosvd_gram  = cell(length(ranks),1);
err_svdsketch  = cell(length(ranks),1);

%% To store train error
err_tr_rhosvd_krp_mttkrps  = zeros(length(ranks),runs);
err_tr_rhosvd_krp  = zeros(length(ranks),runs);
err_tr_rhosvd = zeros(length(ranks),runs);
err_tr_hosvd_gram  = zeros(length(ranks),runs);
err_tr_svdsketch  = zeros(length(ranks),runs);


test_length = size(uctest,4);
rng(1);
% Main loop
for j = 1:length(ranks)
    r = [ranks(j),ranks(j),ranks(j),ranks(j)]; 

    time_rhosvd_krp_memo{j}  = zeros(runs,4);
    time_rhosvd_krp{j} = zeros(runs,4);
    time_rhosvd{j} = zeros(runs,4);
    time_hosvd_gram{j} = zeros(runs,4);
    time_svdsketch{j} = zeros(runs,4);

    err_rhosvd_krp_mttkrps{j} = zeros(runs,test_length);
    err_rhosvd_krp{j} = zeros(runs,test_length);
    err_rhosvd{j} = zeros(runs,test_length);
    err_hosvd_gram{j} = zeros(runs, test_length);
    err_svdsketch{j} = zeros(runs, test_length);


    % SPST-KRP with MTTKRPS (Structure-Preserving HOSVD with KRP and Memoization)
    disp(['RHOSVD-KRP-MEMO (r = ',num2str(ranks(j)),')'])
    for i = 1:runs
        [U,inds,times] = sp_rhosvd_krp_memo(uctrain,r,p,modes,indSelectionMethod);
        time_rhosvd_krp_memo{j}(i,:) = times;

        %train error
        Ts = uctrain(inds{1},inds{2},inds{3},:);
        Tr = ttm(Ts,U,1:3);
        err_tr_rhosvd_krp_mttkrps(j,i) =  norm(uctrain-Tr)/norm(uctrain);


        % test error
        err_sp = zeros(1,test_length);
        for k = 1:test_length
            M = uctest(:,:,:,k)-ucm;
            Ms = M(inds{1},inds{2},inds{3});

            Mr = ttm(Ms,U,1:3);
            err_sp(k) = norm(M-Mr)/norm(M+ucm);
        end
        err_rhosvd_krp_mttkrps{j}(i,:) =  err_sp;
    end


    % SP_HOSVD-KRP (Structure-Preserving HOSVD with KRP)
    disp(['RHOSVD-KRP (r = ',num2str(ranks(j)),')'])
    for i = 1:runs
        [U,inds,times] = sp_rhosvd_krp(uctrain,r,p,modes,indSelectionMethod);
        time_rhosvd_krp{j}(i,:) = times;

         %train error
        Ts = uctrain(inds{1},inds{2},inds{3},:);
        Tr = ttm(Ts,U,1:3);
        err_tr_rhosvd_krp(j,i) =  norm(uctrain-Tr)/norm(uctrain);


        % test error
        err_sp = zeros(1,test_length);
        for k = 1:test_length
            M = uctest(:,:,:,k)-ucm;
            Ms = M(inds{1},inds{2},inds{3});

            Mr = ttm(Ms,U,1:3);
            err_sp(k) = norm(M-Mr)/norm(M+ucm);
        end
        err_rhosvd_krp{j}(i,:) =  err_sp;
    end

    % RHOSVD
    disp(['RHOSVD (r = ',num2str(ranks(j)),')'])
    for i = 1:runs
        [U,inds,times] = sp_rhosvd(uctrain,r,p,modes,indSelectionMethod);
        time_rhosvd{j}(i,:) =  times;

        %train error
        Ts = uctrain(inds{1},inds{2},inds{3},:);
        Tr = ttm(Ts,U,1:3);
        err_tr_rhosvd(j,i) =  norm(uctrain-Tr)/norm(uctrain);


         % test error
        err_rst = zeros(1,test_length);
        for k = 1:test_length
            M = uctest(:,:,:,k)-ucm;
            Ms = M(inds{1},inds{2},inds{3});

            Mr = ttm(Ms,U,1:3);
            err_rst(k) = norm(M-Mr)/norm(M+ucm);
        end
        err_rhosvd{j}(i,:) = err_rst;
    end


    %HOSVD
    disp(['HOSVD-GRAM (r = ',num2str(ranks(j)),')'])
    for i = 1:runs
        [U,inds,times] = sp_hosvd_gram(uctrain,r,modes,indSelectionMethod);
        time_hosvd_gram{j}(i,:) =  times;

        %train error
        Ts = uctrain(inds{1},inds{2},inds{3},:);
        Tr = ttm(Ts,U,1:3);
        err_tr_hosvd_gram(j,i) =  norm(uctrain-Tr)/norm(uctrain);


         % test error
        err_rst = zeros(1,test_length);
        for k = 1:test_length
            M = uctest(:,:,:,k)-ucm;
            Ms = M(inds{1},inds{2},inds{3});

            Mr = ttm(Ms,U,1:3);
            err_rst(k) = norm(M-Mr)/norm(M+ucm);
        end
        err_hosvd_gram{j}(i,:) = err_rst;
    end


    % RHOSVD-SVDSKETCH
    disp(['SVDSKETCH (r = ',num2str(ranks(j)),')'])
    for i = 1:runs
        [U,inds,times] = SVDSKETCH(uctrain,r,p,modes,indSelectionMethod);
        time_svdsketch{j}(i,:) =  times;

        %train error
        Ts = uctrain(inds{1},inds{2},inds{3},:);
        Tr = ttm(Ts,U,1:3);
        err_tr_svdsketch(j,i) =  norm(uctrain-Tr)/norm(uctrain);


         % test error
        err_rst = zeros(1,test_length);
        for k = 1:test_length
            M = uctest(:,:,:,k)-ucm;
            Ms = M(inds{1},inds{2},inds{3});

            Mr = ttm(Ms,U,1:3);
            err_rst(k) = norm(M-Mr)/norm(M+ucm);
        end
        err_svdsketch{j}(i,:) = err_rst;
    end
end
name = "./Experimental_results/tangaroa_"+num2str(indSelectionMethod)+".mat";
save(name,"time_rhosvd_krp_memo","time_rhosvd_krp","time_rhosvd","time_hosvd_gram","time_svdsketch","err_hosvd_gram","err_rhosvd","err_rhosvd_krp","err_rhosvd_krp_mttkrps","err_svdsketch","ranks","indSelectionMethod")


%%  plot 

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

% file_name = "Fig4.png";
% print(gcf, file_name, '-dpng', '-r300')
