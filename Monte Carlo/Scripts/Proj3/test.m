% load('C:\Users\elias\Matlag\MonteCarlo\Monte Carlo\Data\coal_mine_disasters.mat')
load('C:\Users\Seamu\MATLAB\Projects\MonteCarlo\Monte Carlo\Data\coal_mine_disasters.mat')
plot(T(1:751), 1:751)
%%
burninfact = 1.2;
N_real = 10^4;
N = N_real*burninfact; % count in burn in time

d =10;%number of intervals
rho = 0.009;
Phi = 0.5;
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
    if mod(i,N/100)==0
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

start = ceil(N_real*(burninfact-1));
theta = theta(start:end, :);
lambda = lambda(start:end, :);
t = t(start:end, :);

nmbrBlocks = max(N_real/200, 5);

lambdaBlocks = blockify(lambda, nmbrBlocks);
thetaBlocks = blockify(theta, nmbrBlocks);
tBlocks = blockify(t, nmbrBlocks);

%%
figure(1)
for i = 1:d
    plot(1:nmbrBlocks, lambdaBlocks(:, i))
    hold on
end

figure(2)
for i = 1:d+1
    plot(1:nmbrBlocks, tBlocks(:, i))
    hold on
end

figure(3)
plot(1:nmbrBlocks, thetaBlocks)

%%
    a = zeros(10^3, 4);
    for i=1:10^3
        a(i,:) = rwp([1 3 4 7], 0.5);
    end
    mean(a)
    
%%
plot(T(1:751), 1:751)
hold on
mean_lambda = mean(lambda(N_real*(burninfact-1):end,:),1)
mean_t = mean(t(N_real*(burninfact-1):end,:),1)


for i=1:d
    xline(mean_t(i), '--black', {['i=' num2str(i) '     ' 'lambda=' num2str(mean_lambda(i))]});
    
    % Polyfit
end

