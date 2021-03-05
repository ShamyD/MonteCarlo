load('C:\Users\elias\Matlag\MonteCarlo\Monte Carlo\Data\coal_mine_disasters.mat')
plot(T(1:751), 1:751)
%%
burninfact = 1.2;
N_real = 10^6;
N = N_real*burninfact; % count in burn in time

d =10;%number of intervals
rho = 0.002;
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

nmbrBlocks = N_real/200;


%%
    a = zeros(10^3, 4);
    for i=1:10^3
        a(i,:) = rwp([1 3 4 7], 0.5);
    end
    mean(a)
    
%%

mean_lamba = mean(lambda(N_real*(burninfact-1):end,:),1)
mean_t = mean(t(N_real*(burninfact-1):end,:),1)




