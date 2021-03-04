load('C:\Users\Seamu\MATLAB\Projects\MonteCarlo\Monte Carlo\Data\coal_mine_disasters.mat')
plot(T)

%%
burnoutfact = 1.3;
N = 10^1*burnoutfact; % count in burn in time

d = 3;%amount of breakpoints
rho = 0.5;
Phi = 5;

t1 = 1658;
td1 = 1980;
t=zeros(N, d+1);
t(1, :) = linspace(t1, td1, d+1);

lambda = zeros(N, d);
lambda(1,:) = 1;

theta = zeros(N,1);
theta(1) = 1;

for i = 2:N
    tprev = t(i-1, :);
    tnext = MHstep(tprev, rho, tau);
    
    lambdaPrev = lambda(i-1, :);
    thetaPrev = theta(i-1, :);
    [thetaStep, lambdaStep] = Gstep(lambdaPrev, thetaPrev, T, tnext, Phi);
    
    t(i, :) = tnext;
    theta(i,:) = thetaStep;
    lambda(i,:) = lambdaStep; 
    
end

nbrBlocks 


%%
    a = zeros(10^3, 4);
    for i=1:10^3
        a(i,:) = rwp([1 3 4 7], 0.5);
    end
    mean(a)