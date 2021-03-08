function [theta,lambda, t] = hybridSampler(N_real,d, rho, Phi, tau)
    burninfact = 1.2;
    N = N_real*burninfact; % count in burn in time

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
        if mod(i,N/5)==0
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
rej_rate = sum_rej/N % Should be rougly 70%   
   
start = ceil(N_real*(burninfact-1));
theta = theta(start:end, :);
lambda = lambda(start:end, :);
t = t(start:end, :);
 
end