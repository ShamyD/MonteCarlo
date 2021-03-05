

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
    
    wh(i,:) = gumb_inv(r(i,:));
    
    [mu_star(i), beta_star(i)] = est_gumbel(wh(i,:));
    
end

% Skapa delta_mu och delta_beta.
% Skapa left bound and right bound for both.

