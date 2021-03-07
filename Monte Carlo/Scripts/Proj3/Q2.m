

load('C:\Users\elias\Matlag\MonteCarlo\Monte Carlo\Data\atlantic.txt')

%%
[beta, mu] = est_gumbel(atlantic);
B = 200;
n = size(atlantic,1);
r = rand(B,n);
wh = zeros(B,n);
mu_star = zeros(B,1);
beta_star = zeros(B,1);

F_inv_T_star = zeros(B,1);
T = 3*14*100;
F_inv_T = gumb_inv(mu, beta, (1-1/T));


for i=1:B
    
    wh(i,:) = gumb_inv(mu,beta,r(i,:));
    [beta_star(i), mu_star(i)] = est_gumbel(wh(i,:));
    
    % Calc 100-year return value with the estimated parameters.
    F_inv_T_star(i) = gumb_inv(mu_star(i), beta_star(i), (1-1/T));
    
end

mu_diff = sort(mu_star - mu);
beta_diff = sort(beta_star - beta);

F_inv_T_diff = sort(F_inv_T_star - F_inv_T);

alpha = 0.05;

mu_conf = [mu - mu_diff(floor((1-alpha/2)*B)), mu - mu_diff(ceil(alpha/2*B))]
beta_conf = [beta - beta_diff(floor((1-alpha/2)*B)), beta - beta_diff(ceil(alpha/2*B))]

F_inv_T_conf = F_inv_T - F_inv_T_diff(floor((1-alpha)*B))