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
%Probably wrong interval 
N = 10^4;
Fb = wblcdf(25, lambda, k);
Fa = wblcdf(3.5, lambda, k);

for i=1:N
    R = rand(1, 12);
    prob = R.*(Fb-Fa)+Fa;
    V = wblinv(prob, lambda, k);
    power_t(i, :) = P(V);
end
M = N./(Fb-Fa); %Expected number of draws with 0 included
tau_nt = sum(power_t)./M;
tau_nt

sigma_t = std(power_t);
I_t = 2*(1.96./sqrt(M)).*sigma_t

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
tau_ntIS;

sigma_tIS = std(power_t);
I_tIS = 2*(1.96./sqrt(N)).*sigma_tIS

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

%% Prob(P(V)>0)
Fb - Fa


%% E(Ptot)
Ptot = @(x) (1/8)*1.225*pi*(164^2)*(x.^3);

N = 10^4;

for i = 1:N
   pow(i,:) = Ptot(wblrnd(lambda, k));
end

E = mean(pow); %E is meant to be produced analytically

kvot = tau_nis./E; %IS gave the best result
I = I_is/E;
I;

%% Capacity and Availablity

tau = tau_nis;
MAX_P = 9.5e6;
capacity_factor = tau/MAX_P; %REALLY GOOD (especially in winter)

Fb = wblcdf(25, lambda, k);
Fa = wblcdf(3.5, lambda, k);

availability_factor = Fb-Fa %Worse, always lower than 

%% FRÅGOR 
%Skillnad mellan truncerings-algortimerna? - är de samma sak?
%Kan vi trunkera på sättet vi gjorde först?


%IS bättre än AS?
%Rektangelfördelning på inte hela def-mängd?
%Antithetic higher variance than normal monte carlo

%% Del 3

%% Monte Carlo grejs
f = @(x) wblpdf(x, 9.13, 1.96);
F = @(x) wblcdf(x, 9.13, 1.96);
p = 3;
q = 1.5;
alpha = 0.638;
fvec = @(v1, v2) f(v1)*f(v2)*(1+alpha*(1-F(v1)^p)^(q-1)*(1-F(v2)^p)^(q-1)*(F(v1)^p*(1+p*q)-1)*(F(v2)^p*(1+p*q)-1));

x = linspace(0,50, 100);
sur = zeros(length(x), length(x));
for j = 1:length(x)
    for i = 1:length(x)
        sur(i, j) = (P(x(j)) + P(x(i)));%fvec(x(j),x(i));
    end
end

surf(x, x, sur) %peak at 0.012

%%
N=10000;

sur = mvnrnd([12 12], [7 0;0 7], N);

scatter(sur(:,1), sur(:,2))

%%

x = linspace(0,30, 100);

sur = zeros(length(x), length(x));
for j = 1:length(x)
    for i = 1:length(x)
        sur(i, j) = mvnpdf([x(j) x(i)], [12 12], [7 0;0 7]);
    end
end

surf(x, x, sur) %peak at 0.012

%%

x = linspace(0,40, 100);
norm_variance = 7;%7
sur = zeros(length(x), length(x));
for j = 1:length(x)
    for i = 1:length(x)
        sur(i, j) = fvec(x(j),x(i))*(P(x(j)) + P(x(i)))/mvnpdf([x(j) x(i)], [12 12], [norm_variance 0;0 norm_variance]);
    end
end

surf(x, x, sur) %peak at 0.012
%% Uppgift 3A, simulera
N=10000;
points = mvnrnd([12 12], [7 0;0 7], N);

prod = zeros(N,1);

for i=1:length(points)
    prod(i) = fvec(points(i,1),points(i,2))*(P(points(i,1)) + P(points(i,2)))/mvnpdf([points(i,1) points(i,2)], [12 12], [7 0;0 7]);
end


tau_sum = sum(prod)/N; %sums the columns ie the different samples
tau_sum

%sigma = std(power);
%I = 2*(1.96/sqrt(N)).*sigma;

%% Uppgift 3A, simulera rektangel
N=10^4;

prod = zeros(N,1);
points = rand([N,2]);
points = points.*(25-3.5) + 3.5;

for i=1:N
    
    prod(i) = fvec(points(i,1),points(i,2))*(P(points(i,1)) + P(points(i,2)))*(25-3.5)^2;
end


tau_sum = sum(prod)/N; %sums the columns ie the different samples
tau_sum

%sigma = std(power);
%I = 2*(1.96/sqrt(N)).*sigma;

%% Covariance with rejection sampling (AND Variance of sum)
%Bias in estimate due to sampling from only part of set, however the part
%is a CLEAR majority

x = linspace(0,30, 100);

% sur = zeros(length(x), length(x));
% for j = 1:length(x)
%     for i = 1:length(x)
%         sur(i, j) = fvec(x(j),x(i));%/(0.5333*mvnpdf([x(j) x(i)], [12 12], [7 0;0 7]));
%     end
% end

% surf(x, x, sur)

N=10^4;
count = 0;
V = zeros(N,2);
P1 = zeros(N,1);
P2 = zeros(N,1);
Psum = zeros(N,1);
while count<N
  i = count+1;
  X = rand(1,2)*25;
  U = rand();
  if U <= fvec(X(1), X(2))/(0.013) %K and g cancel
      V(count+1, :) = X;
      P1(i) = P(V(i,1));
      P2(i) = P(V(i,2));
      Psum(i) = P1(i)+P2(i);
      count = count + 1
  end
end

C = cov(P1, P2); %variance similar
C
Var = var(Psum)
Std = sqrt(Var)

%% P(Ptot < 9.5)

N = 10^4;

count = 0;
V = zeros(N,2);
fi = zeros(N,1);
Psum = zeros(N,1);
while count<N
  i = count+1;
  X = rand(1,2)*25;
  U = rand();
  if U <= fvec(X(1), X(2))/(0.013) %K and g cancel
      V(count+1, :) = X;
      Psum(i) = P(V(i,1))+P(V(i,2));
      if Psum(i) > 9.5e6
          fi(i) = 1;
      end
      
      count = count + 1
  end
end

prob = mean(fi)
var_fi = var(fi)
Iprob = 2*norminv(0.975)*sqrt(var_fi)/sqrt(N)

