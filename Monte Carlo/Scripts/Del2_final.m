%% Project 1

load powercurve_V164

% SETUP FOR WEIBULL DISTR.
lambda = [10.6, 9.7, 9.2, 8.0, 7.8, 8.1, 7.8, 8.1, 9.1, 9.9, 10.6, 10.6];
k = [2.0, 2.0, 2.0, 1.9, 1.9, 1.9, 1.9, 1.9, 2.0, 1.9, 2.0, 2.0];

%% Standard Monte Carlo
N = 10^4;
power = zeros(N,12);

for i=1:N
    power(i,:) = P(wblrnd(lambda, k));
end

tau_n = sum(power)/N; %sums the columns ie the different samples
tau_n

% 95% Conf. Interval
sigma = std(power);
I = 2*(1.96/sqrt(N)).*sigma;

%% Monte Carlo on truncated distribution (IS framework)
N = 10^4;
Fb = wblcdf(25, lambda, k);
Fa = wblcdf(3.5, lambda, k);

for i=1:N
    R = rand(1, 12);
    prob = R.*(Fb-Fa)+Fa;
    V = wblinv(prob, lambda, k);
    power_t(i, :) = P(V)'.*(Fb-Fa);
end

tau_ntIS = sum(power_t)./N;
tau_ntIS

sigma_tIS = std(power_t);
I_tIS = 2*(1.96./sqrt(N)).*sigma_tIS

%% Importance sampling
N = 10^4;
s = 5; %THIS CAN BE OPTIMISED ()
my = 12;
value = zeros(N,12);

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

%% Antithetic sampling
N = 10^4;
Fb = wblcdf(25, lambda, k);
Fa = wblcdf(3.5, lambda, k);

for i=1:N/2
    U = rand(1, 12);
    Uhat = 1-U;%T(U)
    prob1 = U.*(Fb-Fa)+Fa;
    prob2 = Uhat.*(Fb-Fa)+Fa;
    V1 = P(wblinv(prob1, lambda, k))'.*(Fb-Fa);
    V2 = P(wblinv(prob2, lambda, k))'.*(Fb-Fa);
    
    power(i,:) = (V1 + V2)/2;
end

tau_AS = mean(power);
I_AS = 2*norminv(0.975)*std(power)/sqrt(N/2);

tau_AS
I_AS

tau_n
I

%% Prob(P(V)>0)
Fb - Fa
mean(Fb-Fa)

%% E(Ptot)

EPtot_prim = @(v) (1/8)*1.225*pi*(164^2) * lambda .* exp(-v./k) .* (-k*v - 3*k.^2*v^2 - 6*k.^3*v - 6*k.^4);
EPtot = - EPtot_prim(0)

kvot = tau_nis./EPtot; %IS gave the best result
Ikv = I_is./EPtot;
Ikv
mean(kvot)
%% Capacity and Availablity

tau = tau_nis;
MAX_P = 9.5e6;
capacity_factor = tau/MAX_P; %REALLY GOOD (especially in winter)

Fb = wblcdf(25, lambda, k);
Fa = wblcdf(3.5, lambda, k);

availability_factor = Fb-Fa %Worse, always lower than 

%% Plots
k_average = mean(k);
lambda_average = mean(lambda);

fi = @(x) P(x);
f = @(x) wblpdf(x, lambda_average, k_average);
fff = @(x) fi(x)'.*f(x)./normpdf(x, 5, 12);

x = linspace(0,30, 1000);

figure(1)
plot(x, fi(x))
legend('P(v)', 'Location', 'northwest')
title('Power output as function of wind speed')
xlabel('wind speed v [m/s]')
ylabel('Power output P(v) [W]')

figure(2)
plot(x, fff(x))
legend('\phi(v)f(v)/norm(v)', 'Location', 'northwest')
title('Average manipulated objective function')
xlabel('wind speed v [m/s]')
ylabel('\phi(v)f(v)/norm(v)')

% 
% plot(x, fi(x)'.*f(x))
% hold on
% [M,I] = max(fi(x)'.*f(x));
% my = x(I);
% 
% n = @(y) normpdf(y, my, 5);
% plot(x, M*n(x)/n(my))



%%
l = 10.6;
k = 2;

W = @(x) 1 - exp(-(x/l).^k);
F = @(x) l(1-exp(-x)).^(1/k);

x = linspace(0,40, 100);

plot(x, W(x))
hold on
plot(x, F(x))
