load('C:\Users\elias\Matlag\MonteCarlo\Monte Carlo\Data\coal_mine_disasters.mat')
%load('C:\Users\Seamu\MATLAB\Projects\MonteCarlo\Monte Carlo\Data\coal_mine_disasters.mat')
plot(T(1:751), 1:751)
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
rho = 0.0075
Phi = 3;
N = 10^4;
d = 11;

[theta, lambda, t] = hybridSampler(N, d, rho, Phi, tau);
nmbrBlocks = max(N/100, 5);

lambdaBlocks = blockify(lambda, nmbrBlocks);
thetaBlocks = blockify(theta, nmbrBlocks);
tBlocks = blockify(t, nmbrBlocks);


lambdaVar = var(lambdaBlocks)
thetaVar = var(thetaBlocks)
varT = var(tBlocks)
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
%% Simulate many different d. Rho is adjusted so that an acceptance 
% rate of about 30% is achieved.

rho2 = 0.0072;
rho3 = 0.0053;
rho4 = 0.0035;
rho5 = 0.0053;
rho8 = 0.0065;
rho11 = 0.0075;

Phi = 3;
N = 10^4;
d = [2 3 4 5 8 11];

[theta2, lambda2, t2] = hybridSampler(N, 2, rho2, Phi, tau);
[theta3, lambda3, t3] = hybridSampler(N, 3, rho3, Phi, tau);
[theta4, lambda4, t4] = hybridSampler(N, 4, rho4, Phi, tau);
[theta5, lambda5, t5] = hybridSampler(N, 5, rho5, Phi, tau);
[theta8, lambda8, t8] = hybridSampler(N, 8, rho8, Phi, tau);
[theta11, lambda11, t11] = hybridSampler(N, 11, rho11, Phi, tau);
%%
nmbrBlocks = max(N/100, 5);

figure(1)
thetaBlocks2 = blockify(theta2, nmbrBlocks);
plot(thetaBlocks2)
figure(2)
thetaBlocks3 = blockify(theta3, nmbrBlocks);
plot(thetaBlocks3)
figure(3)
thetaBlocks4 = blockify(theta4, nmbrBlocks);
plot(thetaBlocks4)
figure(4)
thetaBlocks5 = blockify(theta5, nmbrBlocks);
plot(thetaBlocks5)
figure(5)
thetaBlocks8 = blockify(theta8, nmbrBlocks);
plot(thetaBlocks8)
figure(6)
thetaBlocks11 = blockify(theta11, nmbrBlocks);
plot(thetaBlocks11)

%% Check theta blck variance dependending on some ds.

hold on
thetaBlocks2 = blockify(theta2, nmbrBlocks);
plot(thetaBlocks2)
thetaBlocks4 = blockify(theta4, nmbrBlocks);
plot(thetaBlocks4)
thetaBlocks11 = blockify(theta11, nmbrBlocks);
plot(thetaBlocks11)

xlabel('Block Sample')
ylabel('Block mean of \theta')
title('Block mean of \theta for d=2,4,11')
legend('d=2', 'd=4', 'd=11', 'Location', 'northeast')

%% Mean and (block)variance for theta for d=2,3,4,5,8,11.
means_theta = [mean(thetaBlocks2) mean(thetaBlocks3) mean(thetaBlocks4) mean(thetaBlocks5) mean(thetaBlocks8) mean(thetaBlocks11)]
var_theta = [var(thetaBlocks2) var(thetaBlocks3) var(thetaBlocks4) var(thetaBlocks5) var(thetaBlocks8) var(thetaBlocks11)]

%% Plot t over samples for d=8.
tBlocks8 = blockify(t8, nmbrBlocks);
plot(tBlocks8)
var(tBlocks8)

xlabel('Block Sample')
ylabel('Block mean of breakpoints t_i')
title('Block mean of t_i for d=8')
legend('i=1', 'i=2', 'i=3', 'i=4', 'i=5', 'i=6', 'i=7', 'i=8', 'i=9', 'Location', 'southeast')

%%

for subf=1:6
    figure(subf)
    plot(T(1:751), 1:751)
    hold on
    if subf==1
        mean_lambda = mean(lambda2,1)
        mean_t = mean(t2,1)
        xlabel('Time (years)')
        ylabel('Cumulative distribution of accidents')
        title(['Breakpoints for d=' num2str(d(subf))])
    end
    if subf==2
        mean_lambda = mean(lambda3,1)
        mean_t = mean(t3,1)
        xlabel('Time (years)')
        ylabel('Cumulative distribution of accidents')
        title(['Breakpoints for d=' num2str(d(subf))])
    end
    if subf==3
        mean_lambda = mean(lambda4,1)
        mean_t = mean(t4,1)
        xlabel('Time (years)')
        ylabel('Cumulative distribution of accidents')
        title(['Breakpoints for d=' num2str(d(subf))])
    end
    if subf==4
        mean_lambda = mean(lambda5,1)
        mean_t = mean(t5,1)
        xlabel('Time (years)')
        ylabel('Cumulative distribution of accidents')
        title(['Breakpoints for d=' num2str(d(subf))])
    end
    if subf==5
        mean_lambda = mean(lambda8,1)
        mean_t = mean(t8,1)
        xlabel('Time (years)')
        ylabel('Cumulative distribution of accidents')
        title(['Breakpoints for d=' num2str(d(subf))])
    end
    if subf==6
        mean_lambda = mean(lambda11,1)
        mean_t = mean(t11,1)
        xlabel('Time (years)')
        ylabel('Cumulative distribution of accidents')
        title(['Breakpoints for d=' num2str(d(subf))])
    end


    for i=1:d(subf)
        xline(mean_t(i), '--black', {['i=' num2str(i) '     ' 'lambda=' num2str(mean_lambda(i))]});

        % Polyfit
    end
end

% - Simulate from our best approximation of the parameters. See if data
% similar.
% - Make grid of time and generate mean of lambdas for each gridpoint.


%%
plot(T(1:751), 1:751)
hold on
mean_lambda = mean(lambda,1)
mean_t = mean(t,1)


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
lambda_for_grid = zeros(N,inner_grid_points);
% 
% lambda_wo_burnin = lambda(N-N_real:end,:);
% t_wo_burnin = t(N-N_real:end,:);

for i=1:N
    lambda_for_grid(i,:) = lambda_for_tp(lambda(i,:),d,t(i,:), grid);
end

lambda_mean = mean(lambda_for_grid,1);

plot(grid(2:end-1), lambda_mean)



%%
% Plot using time-grid for MANY DIFFERENT d.
for subf=1:6    
    grid_size = 0.1;
    grid = 1658:grid_size:1980;
    inner_grid_points = (1980-1658)/grid_size - 1;
    lambda_for_grid = zeros(N,inner_grid_points);
    
    if subf==1
        lambda_d = lambda2;
        t_d = t2;
    end
    if subf==2
        lambda_d = lambda3;
        t_d = t3;
    end
    if subf==3
        lambda_d = lambda4;
        t_d = t4;
    end
    if subf==4
        lambda_d = lambda5;
        t_d = t5;
    end
    if subf==5
        lambda_d = lambda8;
        t_d = t8;
    end
    if subf==6
        lambda_d = lambda11;
        t_d = t11;
    end
    
    
    for i=1:N
        lambda_for_grid(i,:) = lambda_for_tp(lambda_d(i,:),d(subf),t_d(i,:), grid);
    end

    lambda_mean = mean(lambda_for_grid,1);
    
    figure(subf)
    plot(grid(2:end-1), lambda_mean)
    xlabel('Time (years)')
    ylabel('Disaster Intensity (1/year)')
    title(['Average Intensity for d=' num2str(d(subf))])

end






