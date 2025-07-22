clc
clear all
close all

maxNumCompThreads(1);

addpath("sysID/")
load mat/data_155_oct19.mat

%setting parameter values
ts = .01;
ni = size(B,2);
no = size(C,1);

p = 200;
q = p;
N = 2*p-1;

runs = 10;
Ranks = [75,100,125,150];


Time_rsvd = zeros(length(Ranks),runs);
Time_gauss = zeros(length(Ranks),runs);
Time_krp = zeros(length(Ranks),runs);

hd_rsvd = zeros(length(Ranks),runs);
hd_gauss = zeros(length(Ranks),runs);
hd_krp = zeros(length(Ranks),runs);



%% convert to discrete
sys = ss(A,B,C,0);
sysd = c2d(sys,ts,'tustin');
[Ad,Bd,Cd,Dd] = ssdata(sysd);

eig_Ad = eig(Ad);


% Markov parameters
markov = cell(1,N);
T = zeros(no,N,ni);
f = Bd;
fact = [1:p-1,p:-1:1];
for jj = 1:N
    g = Cd*f;
    f = Ad*f;           
    markov{jj} = g;
    T(:,jj,:) = sqrt(jj)*g;
end

% Setup matrices
mat.E = cell(1,N);
mat.M = markov;

I = eye(p,p);
for j = 1:p
    mat.E{j} = sparse(hankel(I(:,j)));
end
for j = 2:p
    mat.E{p+j-1} = sparse(hankel(zeros(p,1), I(j,:)));
end

%% Compare blockstructure
x1 = randn(q,1);
x2 = randn(ni,1);
x = kr(x1,x2);
Hf = makehankelfun(markov,p,q,no,ni,'notransp');
y1 = Hf(x, 'notransp');

y2 = 0*y1;
for j = 1:2*p-1
        y2= y2  + kr(mat.E{j}*x1,mat.M{j}*x2);
end
assert(norm(y1-y2)/norm(y1)<1e-14)


%% Running time experiments
for k =1:length(Ranks)
    r = Ranks(k);
    disp(['Single-View Randomized SVD with Dense Gaussian at r = ',num2str(r)])
    for i = 1:runs
        t1 = tic;
        [U1,S1,V1] = blockstructuredsvd_dense(mat, r);
        [A1,B1,C1] = era(U1,S1,V1,r,ni,no);
        Time_gauss(k,i) = toc(t1);
        lr1 = eig(A1);
        [hd,~] = HausdorffDist(lr1,eig_Ad);
        hd_gauss(k,i) = hd;
    end
    
    disp(['Single-View Randomized SVD with KRP structured random matrices at r = ',num2str(r)])
    for i = 1:runs
        t2 = tic;
        [U2,S2,V2] = blockstructuredsvd_krp(mat, r);
        [A2,B2,C2] = era(U2,S2,V2,r,ni,no);
        Time_krp(k,i) = toc(t2);
        lr2 = eig(A2);
        [hd,~] = HausdorffDist(lr2,eig_Ad);
        hd_krp(k,i) = hd;
    end
        
    disp(['RandSVD at r = ',num2str(r)])
    for i = 1:runs
        t3 = tic;
        [A3,B3,C3,~] = impulse_era(markov,p,q,no,ni,r,'randsvd'); 
        Time_rsvd(k,i) = toc(t3);
        lr3 = eig(A3);
        [hd,~] = HausdorffDist(lr3,eig_Ad);
        hd_rsvd(k,i) = hd;
    end
end

file_name = './Experimental_results/SysID_results.mat';
save(file_name, 'Time_rsvd', 'Time_gauss', 'Time_krp', 'hd_rsvd', 'hd_gauss', 'hd_krp', 'Ranks');


%% Ploting
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
set(gcf, 'Position', [100, 100, 900, 300]);  % [left, bottom, width, height]
% ---- Subplot 1: Hausdorff Distance ----
subplot(1,2,2)
plot(Ranks, abs(avg_hd_krp),   'r*-', 'DisplayName', 'R-KRP',    'LineWidth', 2); hold on 
plot(Ranks, abs(avg_hd_gauss), 'ks-', 'DisplayName', 'R-Dense',  'LineWidth', 2); 
plot(Ranks, abs(avg_hd_rsvd),   'bo-', 'DisplayName', 'RandSVD',     'LineWidth', 2); 

set(gca, 'FontSize', 14.8);
xlabel('Rank (r)');
ylabel('Average Hausdorff Distance');
legend('show', 'Location', 'northeast');
title('Error')
xticks(Ranks);


% ---- Subplot 2: Average Runtime ----
subplot(1,2,1)
plot(Ranks, avg_krp,    'r*-', 'DisplayName', 'R-KRP',    'LineWidth', 2); hold on
plot(Ranks, avg_gauss,  'ks-', 'DisplayName', 'R-Dense',  'LineWidth', 2); 
plot(Ranks, avg_rsvd,   'bo-', 'DisplayName', 'RandSVD',     'LineWidth', 2); 

set(gca, 'FontSize', 14.8);
xlabel('Rank (r)');
ylabel('Average Runtime (Seconds)');
legend('show', 'Location', 'northwest');
xticks(Ranks);
title('Runtime')

% file_name = "Fig1.png";
% print(gcf, file_name, '-dpng', '-r300');

%% Number of random numbers 
for i = 1:length(Ranks)
    r = Ranks(i);
    [m,n] = size(mat.M{1});
    [p,q] = size(mat.E{1});
    rng1 = (n*q)*(r+20);
    rng2 = (n+q)*(2*r+1) + (m+p)*(2*(2*r+1)+1);
    rng3 = (n*q)*(2*r+1) + (m*p)*(2*(2*r+1)+1);
    % rng2/rng1;
    disp([num2str(r), ' ', num2str((rng2/rng1)*100), ' ', num2str((rng2/rng3)*100)]);
end

