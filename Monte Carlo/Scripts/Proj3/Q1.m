% load('C:\Users\elias\Matlag\MonteCarlo\Monte Carlo\Data\coal_mine_disasters.mat')
load('C:\Users\Seamu\MATLAB\Projects\MonteCarlo\Monte Carlo\Data\coal_mine_disasters.mat')
plot(T(1:751), 1:751)
xlabel('year')
ylabel('Number of disasters since 1658')
title('Coal mine disasters between 1658 and 1980') 
tau = T';
%% Hybrid sampler - now replaced with function
burninfact = 1.2;
N_real = 10^4;
N = N_real*burninfact; % count in burn in time

d =10;%number of intervals
rho = 0.009;
Phi = 1;
tau = T';

t1 = 1658;
td1 = 1980;
t=zeros(N, d+1);
t(1, :) = linspace(t1, td1, d+1);
ni = create_ni(tau, t(1,:));

lambda = zeros(N, d);
lambda(1,:) = 1;

theta = zeros(N,1);
theta(1) = 1;

sum_rej = 0;
for i = 2:N
    if mod(i,N/1000)==0
        i
    end
    tprev = t(i-1, :);
    lambdaPrev = lambda(i-1, :);
    
    [tnext, ni, sum_rej] = MHstep(tprev, rho, tau, lambdaPrev, ni, sum_rej, i);
    [thetaStep, lambdaStep] = Gstep(lambdaPrev, ni, tnext, Phi, d);
    
    t(i, :) = tnext;
    theta(i,:) = thetaStep;
    lambda(i,:) = lambdaStep; 
    
end

%% Produce blocks and remove burn in
rho = 0.0055;
Phi = 200;
N = 10^5;
d = 4;

[theta, lambda, t] = hybridSampler(N, d, rho, Phi, tau);
nmbrBlocks = max(N/1000, 5);

lambdaBlocks = blockify(lambda, nmbrBlocks);
thetaBlocks = blockify(theta, nmbrBlocks);
tBlocks = blockify(t, nmbrBlocks);

lambdaVar = var(lambdaBlocks)
lambdaMean = mean(lambdaBlocks,1)
thetaVar = var(thetaBlocks)
thetaMean = mean(thetaBlocks, 1)
tVar = var(tBlocks)
tMean = mean(tBlocks,1)
%%
[~,d] = size(lambdaBlocks);

lambdaLegend = cell([1 d]);
for i = 1:d
    lambdaLegend{i} = append('\lambda', int2str(i));
end

figure(1)
plot(1:nmbrBlocks, lambdaBlocks)
xlabel('Block b')
ylabel('Intensity \lambda')
title('Intensities \lambda for different blocks')
legend(lambdaLegend, 'Location', 'northwest')

tLegend = cell([1 d+1]);
for i = 1:d+1
    tLegend{i} = append('t', int2str(i));
end

figure(2)
plot(1:nmbrBlocks, tBlocks)
xlabel('Block b')
ylabel('breakpoint t')
title('Breakpoints t for different blocks')
legend(tLegend, 'Location', 'northwest')

figure(3)
plot(1:nmbrBlocks, thetaBlocks)
xlabel('Block b')
ylabel('Parameter \theta')
title('Parameter \theta for different blocks')    

figure(4)
plot(autocorr(tBlocks(:,2)))
hold on
plot(autocorr(t(:,2)))
legend('block', 'nonBlock')
%%
plot(T(1:751), 1:751)
hold on
mean_lambda = mean(lambda(N_real*(burninfact-1):end,:),1)
mean_t = mean(t(N_real*(burninfact-1):end,:),1)


for i=1:d
    xline(mean_t(i), '--black', {['i=' num2str(i) '     ' 'lambda=' num2str(mean_lambda(i))]});
    
    % Polyfit
end

% - Simulate from our best approximation of the parameters. See if data
% similar.
% - Make grid of time and generate mean of lambdas for each gridpoint.

%%
% Plot using time-grid.
size = 1;
grid = 1658:size:1980;
inner_grid_points = (1980-1658)/size - 1;
lambda_for_grid = zeros(N_real,inner_grid_points);

lambda_wo_burnin = lambda(N-N_real:end,:);
t_wo_burnin = t(N-N_real:end,:);

for i=1:N_real
    lambda_for_grid(i,:) = lambda_for_tp(lambda_wo_burnin(i,:),d,t_wo_burnin(i,:), grid);
end

lambda_mean = mean(lambda_for_grid,1);

plot(grid(2:end-1), lambda_mean)



%% Produce plots for different rho:s

rhoBold = 0.07;
rhoReg = 0.0055;
rhoCaut = 0.0005;
Phi = 3;
N = 10^5;
d = 4;

nmbrBlocks = max(N/1000, 5);

[theta, lambda, t] = hybridSampler(N, d, rhoBold, Phi, tau);

lambdaBlocksBold = blockify(lambda, nmbrBlocks);
thetaBlocksBold = blockify(theta, nmbrBlocks);
tBlocksBold = blockify(t, nmbrBlocks);

[theta, lambda, t] = hybridSampler(N, d, rhoReg, Phi, tau);

lambdaBlocksReg = blockify(lambda, nmbrBlocks);
thetaBlocksReg = blockify(theta, nmbrBlocks);
tBlocksReg = blockify(t, nmbrBlocks);

[theta, lambda, t] = hybridSampler(N, d, rhoCaut, Phi, tau);

lambdaBlocksCaut = blockify(lambda, nmbrBlocks);
thetaBlocksCaut = blockify(theta, nmbrBlocks);
tBlocksCaut = blockify(t, nmbrBlocks);

%Results for each rho
lambdaVarBold = var(lambdaBlocksBold);
lambdaMeanBold = mean(lambdaBlocksBold,1);
thetaVarBold = var(thetaBlocksBold);
thetaMeanBold = mean(thetaBlocksBold, 1);
tVarBold = var(tBlocksBold);
tMeanBold = mean(tBlocksBold,1);

lambdaVarReg = var(lambdaBlocksReg);
lambdaMeanReg = mean(lambdaBlocksReg,1);
thetaVarReg = var(thetaBlocksReg);
thetaMeanReg = mean(thetaBlocksReg, 1);
tVarReg = var(tBlocksReg);
tMeanReg = mean(tBlocksReg,1);

lambdaVarCaut = var(lambdaBlocksCaut);
lambdaMeanCaut = mean(lambdaBlocksCaut,1);
thetaVarCaut = var(thetaBlocksCaut);
thetaMeanCaut = mean(thetaBlocksCaut, 1);
tVarCaut = var(tBlocksCaut);
tMeanCaut = mean(tBlocksCaut,1);

%%
[~,d] = size(lambdaBlocks);


figure(3)
plot(1:nmbrBlocks, thetaBlocks)
xlabel('Block b')
ylabel('Parameter \theta')
title('Parameter \theta for different blocks')    

%%

plot(autocorr(tBlocksBold(:,3)))
hold on
plot(autocorr(tBlocksReg(:,3)))
hold on
plot(autocorr(tBlocksCaut(:,3)))
legend('Bold', 'Regular', 'Cautious')
xlabel('Block difference (l)')
ylabel('Auttocorrelation \rho(l)')
title('Autocorrelation of t3 for Block differences l ')    

%%

tLegend = cell([1 d+1]);
for i = 1:d+1
    tLegend{i} = append('t', int2str(i));
end

figure(2)
plot(1:nmbrBlocks, tBlocksCaut)
xlabel('Block b')
ylabel('breakpoint t')
title('Breakpoints t for different blocks, for cautious choice of \rho')
legend(tLegend, 'Location', 'northwest')

%%
lambdaLegend = cell([1 d]);
for i = 1:d
    lambdaLegend{i} = append('\lambda', int2str(i));
end

figure(1)
plot(1:nmbrBlocks, lambdaBlocksCaut)
xlabel('Block b')
ylabel('Intensity \lambda')
title('Intensities \lambda for different blocks with cautious choice of \rho')
legend(lambdaLegend, 'Location', 'northwest')
