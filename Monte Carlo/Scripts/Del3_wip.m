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
%% Uppgift 3A, simulera - IRRELEVANT
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

%% Uppgift 3A, simulera rektangel - IRRELEVANT
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

%% Covariance with rejection sampling (AND Variance of sum) - IRRELEVANT
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

%% P(Ptot < 9.5) - IRRELEVANT

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

%% 3a
%a E(P(v1) + P(v2)) = E(P(v1)) + E(P(v2)) = tau_nis
%Behlver nya k och lambda

lambda = 9.13;
k = 1.96;

N = 10^4;
s = 5; %THIS CAN BE OPTIMISED ()
my = 12;
value = zeros(N,1);

for i=1:N
   X = normrnd(my, s, 1);
   value(i) = wblpdf(X, lambda, k).*P(X)'./normpdf(X, my, s);
end

tau_nis2 = sum(value)/N;
tau_sum = 2*tau_nis2;   


%% 3b-c
f = @(x) wblpdf(x, lambda, k);
a = 3.5;
b = 25;

N = 10^4;
Pprod = zeros(N,1);
C = zeros(N,1);
V = rand(N,2)*(b-a) +a;
P2 = zeros(N,1);

for i=1:N
    Pprod(i) = P(V(i,1))*P(V(i,2))*f(V(i,1))*fvec(V(i,1),V(i,2))*(b-a)^2;
    P2(i) = P(V(i,1))^2*f(V(i,1))*(b-a); %calc 3c
    C(i) = f(V(i,1))*fvec(V(i,1), V(i,2))*(b-a)^2; %Probably not nececarry
end %Should try not to divide up

P1P2 = mean(Pprod)/mean(C);

covariance = P1P2 - tau_nis2^2;

Ep2 = mean(P2);
varp = (Ep2 - tau_nis2^2);
varsum = 2*varp + 2*covariance;%result 3c
stdsum = sqrt(varsum);%result 3c

%% Plot 3d

x = linspace(25, 100, 100);
func = zeros(length(x), 1);
v2 = 14;
for i = 1:length(x)
func(i) = fvec(x(i), v2)/(2*normpdf(x(i), 25, 6));
end
plot(x, func)

%% Psum = 9.5
N = 10^4;
v1 = rand(N, 1)*(25-14)+14;
m = 14;
sig = 6;
v2 = normrnd(m, sig, N, 1);
for i =1:N %flip normal dist
    v = v2(i);
    if v < m
        v2(i) = 2*m - v;
    end
end

p95 = zeros(N,1);
for i = 1:N
    p95 = fvec(v1(i), v2(i))*(25-14)/(2*normpdf(v2(i), m, sig));
end

prob95 = mean(p95);
tot_prob95 = 2*prob95

%% 3d
N = 10^4;
a = 3.5;
b = 25;
V = rand(N, 2)*(b-a)+a;
Pg95 = zeros(N,1);

for i = 1:N
    psum = P(V(i,1)) + P(V(i,2));
    if psum > 9.5
        Pg95(i) = fvec(V(i,1), V(i,2))*(b-a)^2;
    end
end

prob95g = mean(Pg95)
prob95l = 1 - prob95g - tot_prob95

%% 3d sampling from R2+
N = 10^4;
g = @(x) wblpdf(x, lambda, k);
V = wblrnd(lambda, k, N, 2);

Psumfi = zeros(N,2); %k1:greater than

for i = 1:N
    if P(V(i, 1)) + P(V(i,2)) > 9.5e6
        Psumfi(i, 1) = fvec(V(i,1), V(i,2))/(g(V(i,1))*g(V(i,2)));
    elseif P(V(i, 1)) + P(V(i,2)) < 9.5e6
        Psumfi(i, 2) = fvec(V(i,1), V(i,2))/(g(V(i,1))*g(V(i,2)));
    end
end

pgreat = mean(Psumfi(:,1))
pless = mean(Psumfi(:,2))
diff = 1 - pgreat - pless 

Iplus = 2*norminv(0.975)*sqrt(pgreat*(1-pgreat)/N)
Iminus = 2*norminv(0.975)*sqrt(pless*(1-pless)/N)