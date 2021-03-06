

load('C:\Users\elias\Matlag\MonteCarlo\Monte Carlo\Data\atlantic.txt')

%%

[beta, mu] = est_gumbel(atlantic);
B = 200;
n = size(atlantic,1);
r = rand(B,n);
wh = zeros(B,n);
mu_star = zeros(B,1);
beta_star = zeros(B,1);

for i=1:B
    
    wh(i,:) = gumb_inv(mu,beta,r(i,:));
    [beta_star(i), mu_star(i)] = est_gumbel(wh(i,:));
    
    % Calc 100-year return value with the estimated parameters.
    F_inv_T_star(i) = gumb_inv(mu_star(i), beta_star(i), (1-1/T));
    
end

mu_diff = sort(mu_star - mu);
beta_diff = sort(beta_star - beta);

