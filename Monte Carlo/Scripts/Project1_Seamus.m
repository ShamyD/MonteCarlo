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

%% Monte Carlo on Truncated distribution
N = 10^4;
Fb = wblcdf(25, lambda, k);
Fa = wblcdf(3.5, lambda, k);

for i=1:N
    R = rand(1, 12);
    prob = R.*(Fb-Fa)+Fa;
    V = wblinv(prob, lambda, k);
    power_t(i, :) = P(V);
end
M = N./(Fb-Fa) %Expected number of draws with 0 included
tau_nt = sum(power_t)./M;
tau_nt

sigma_t = std(power_t);
I_t = 2*(1.96./sqrt(M)).*sigma_t;

%% Importance sampling
N = 10^4;
s = 5; %THIS CAN BE OPTIMISED ()
my = 12;

for i=1:N
   X = normrnd(my, s, 1, 12);
   value(i,:) = wblpdf(X, lambda, k).*P(X)'./normpdf(X, my, s);
end

tau_nis = sum(value)/N;
tau_nis

sigma_is = std(value);
I_is = 2*(1.96./sqrt(N)).*sigma_is;
I_is

% SUGGESTION
% TRUNCATED (a = 3.5, b = 25) NORMAL DISTRIBUTION WITH MEAN IN EXPECTED
% VALUE FOR EVERY MONTH

%% Försök till att hitta optimal
fi = @(x) P(x);
f = @(x) wblpdf(x, 10.6, 2.0);
fff = @(x) fi(x)'.*f(x)+1;

x = linspace(0,30, 1000);
plot(x, fi(x)'.*f(x))
hold on
[M,I] = max(fi(x)'.*f(x));
my = x(I);

n = @(y) normpdf(y, my, 5);
% plot(x, M*n(x)/n(my))
%%
plot(x, fi(x)'.*f(x)./n(x))



