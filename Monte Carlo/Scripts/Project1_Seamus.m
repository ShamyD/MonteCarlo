%% Project 1

load powercurve_V164

%% Test P-function
P([6.5, 5.6])

%% SETUP FOR WEIBULL DISTR.
lambda = [10.6, 9.7, 9.2, 8.0, 7.8, 8.1, 7.8, 8.1, 9.1, 9.9, 10.6, 10.6];
k = [2.0, 2.0, 2.0, 1.9, 1.9, 1.9, 1.9, 1.9, 2.0, 1.9, 2.0, 2.0];

V = wblrnd(lambda, k)
%% Standard Monte Carlo
N = 10^4;

for i=1:N
    power(i,:) = P(wblrnd(lambda, k));
end

tau_n = sum(power)/N; %sums the columns ie the different samples
tau_n

%% 95% Conf. Interval
sigma = std(power);
I = 2*(1.96/sqrt(N)).*sigma;