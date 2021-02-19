%% Part 3 - Setup - Run sections in order

load powercurve_V164
lambda = 9.13;
k = 1.96;
f = @(x) wblpdf(x, 9.13, 1.96);
F = @(x) wblcdf(x, 9.13, 1.96);
p = 3;
q = 1.5;
alpha = 0.638;
fvec = @(v1, v2) f(v1)*f(v2)*(1+alpha*(1-F(v1)^p)^(q-1)*(1-F(v2)^p)^(q-1)*(F(v1)^p*(1+p*q)-1)*(F(v2)^p*(1+p*q)-1));


%% 3a
%a E(P(v1) + P(v2)) = E(P(v1)) + E(P(v2)) = tau_nis

N = 10^4;
s = 5; %Can be improved
my = 12;
value = zeros(N,1);

for i=1:N
   X = normrnd(my, s, 1);
   value(i) = wblpdf(X, lambda, k).*P(X)'./normpdf(X, my, s);
end %IMPORTANCE SAMPLING SIMILAR TO 2b

tau_nis2 = sum(value)/N; %Expectation of P(v1)
tau_sum = 2*tau_nis2; %Expectation of the sum


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
    C(i) = f(V(i,1))*fvec(V(i,1), V(i,2))*(b-a)^2; %For self normalized importance sampling
end %Should try not to divide up

P1P2 = mean(Pprod)/mean(C); %E(P(v1)*P(v2))
covariance = P1P2 - tau_nis2^2; %C(P(v1),P(v2))

Ep2 = mean(P2); %E(p(v1)^2)
varp = (Ep2 - tau_nis2^2); %E(p(v1)^2) - E(p(v1))^2 = V(P(v1))
varsum = 2*varp + 2*covariance;% V(P(v1)+P(v2))
stdsum = sqrt(varsum);%D(P(v1)+P(v2))

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

pgreat = mean(Psumfi(:,1)) %P(Psum > 9.5)
pless = mean(Psumfi(:,2)) %P(Psum < 9.5)
diff = 1 - pgreat - pless 

Iplus = [pgreat - norminv(0.975)*sqrt(pgreat*(1-pgreat)/N) pgreat + norminv(0.975)*sqrt(pgreat*(1-pgreat)/N)] %Interval for pgreat
Iminus = [pless - norminv(0.975)*sqrt(pless*(1-pless)/N) pless + norminv(0.975)*sqrt(pless*(1-pless)/N)]

%% Plots - plots that were used when studying functions

figure(1)
x = linspace(25, 100, 100);
func = zeros(length(x), 1);
v2 = 14;
for i = 1:length(x)
func(i) = fvec(x(i), v2)/(2*normpdf(x(i), 25, 6));
end
plot(x, func)

figure(2)
x = linspace(0,40, 100);
norm_variance = 7;%7
sur = zeros(length(x), length(x));
for j = 1:length(x)
    for i = 1:length(x)
        sur(i, j) = fvec(x(j),x(i))*(P(x(j)) + P(x(i)))/mvnpdf([x(j) x(i)], [12 12], [norm_variance 0;0 norm_variance]);
    end
end
surf(x, x, sur) %peak at 0.012

figure(3)
x = linspace(0,50, 100);
sur = zeros(length(x), length(x));
for j = 1:length(x)
    for i = 1:length(x)
        sur(i, j) = (P(x(j)) + P(x(i)));%fvec(x(j),x(i));
    end
end
surf(x, x, sur) %peak at 0.012

