


% Q 3.

% Generate N walks from g(x_(0:n))
N=1000;
n=20;
d = 2;
walks = zeros(N,d*(n + 1));
sa = zeros(N,n+1);

for c = 1:n
    % Generate new points.
    prev_col = walks(:,2*c - 1:2*c);
    new_col = getNewStep(N,d) + prev_col;
    walks(:,d*c + 1:d*c + d) = new_col;
    
    %Loop over previous points to see which have already been
    %visited.
    memTot = zeros(N,c);
    for i=1:c
        old_col = walks(:, d*i - d + 1:d*i);
        temp_col = old_col - new_col;
        members = ismember(temp_col,zeros(1,d), 'rows');
        memTot(:,i) = members;
    end    
    sa(:,c+1) = sum(memTot')';
end


sa_count = zeros(N,n+1);
sa_count(:,3) = sa(:,3);
for j=4:n+1
    sa_count(:,j) = sa_count(:,j-1) + sa(:,j);
end

% Sum up number of self-avoiding walks Nsa for n=1:n.
Nsa = sum(sa_count==0);

% Estimate c_n:s from Nsa.
fracerino = Nsa / N;
Narray = 4.^(linspace(0,n,n+1));
cns = fracerino.*Narray


%Plot.
semilogy(0:n,cns)
xlabel('Steps (n)')
ylabel('Nr of self-avoiding walks, c_n')
title('Naive approach to simulating RW:s')

%% Q 4


% Generate N walks from g(x_(0:n))
N=1000;
n = 20;
d = 2;
walks = zeros(N,d*(n + 1));
N_tot = ones(N,n);
G = 4*ones(N,n);
result = ones(1,n);
dirs = getDirs(d);

for c = 1:n
    
    
    
    % Get which free nb:s there are and how many (for all walks).
    [free_nb, nr_free_nb] = getFreeNb(walks(:,1:d*c),N,d,dirs);
    
    % Get new point for each of the N separate walks.
    prev_col = walks(:,d*c - d + 1:d*c);
    new_col = getFreeStep(N,d,prev_col, free_nb, nr_free_nb);
    walks(:,d*c + 1:d*c + d) = new_col;
    
    %Update nr_free_nb_tot
%     N_tot(:, c) = nr_free_nb;
%     invgn = prod(N_tot(:,1:c)', 1);
%     result(c) = sum(invgn)/N;
%

    N_tot(:, c) = nr_free_nb;
    if c>1
        G(:,c) = nr_free_nb.*G(:,c-1);
    end
end

result = sum(G,1)/N

%Plot.
semilogy(1:n,result)
xlabel('Steps (n)')
ylabel('Nr of self-avoiding walks, c_n')
title('Only self-avoiding RW:s generated')

%% Q5 and onwards
N = 1000;
n = 20;
d = 2;
rep = 1;
cum_sums = zeros(rep, n);

for i = 1:rep
    cum_sums(i, :) = SISR_sample(N, n, d);
    if mod(i, 5) == 0
        i
    end
end

c_SISR_mean = mean(cum_sums,1)


%Plot.
semilogy(1:n,c_SISR_mean)
xlabel('Steps (n)')
ylabel('Nr of self-avoiding walks, c_n')
title('SISR approach')

% figure(2);
% plot(1:n, c_SISR_mean)
%%
getParams(cum_sums);

%% Final questions
N = 1000;
n = 30;
rep = 30;
dimensions = 10;

cum_sums = zeros(rep, n);
C_means = zeros(dimensions, n);

mu_means = zeros(dimensions,1);
mu_bounds = zeros(dimensions, rep);

gamma_means = zeros(dimensions,1);

A_means = zeros(dimensions,1);
A_bounds = zeros(dimensions, rep);

%Loop over all dimensions tested
for d = 1:dimensions
    cum_sums = SISR_sampling(N, n, d, rep);
    C_means(d, :) = mean(cum_sums,1);
    [mu, gamma, A] = getParams(cum_sums);
    
    mu_means(d,1) = mean(mu,1);
    mu_bounds(d, :) = checkBounds(mu,d, 2*d-1)';
    
    gamma_means(d,1) = mean(gamma,1);
    
    A_means(d,1) = mean(A,1);
    A_bounds(d,:) = A>=1';
    
end


