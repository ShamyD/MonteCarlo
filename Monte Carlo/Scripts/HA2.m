


% Q 3.

% Generate N walks from g(x_(0:n))
N=10000;
n=10;
d = 2;
walks = zeros(N,d*(n + 1));
sa = zeros(N,n+1);

for c = 1:n
    
    prev_col = walks(:,2*c - 1:2*c);
    new_col = getNewStep(N,d) + prev_col;
    walks(:,d*c + 1:d*c + d) = new_col;
    
    
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

Nsa = sum(sa_count==0);

fracerino = Nsa / N;
Narray = 4.^(linspace(0,n,n+1));
cns = fracerino.*Narray

%% Q 4


% Generate N walks from g(x_(0:n))
N=10000;
n = 10;
d = 2;
walks = zeros(N,d*(n + 1));
N_tot = ones(N,n);
G = 4*ones(N,n);
result = ones(1,n);

for c = 1:n
    % Get which free nb:s there are and how many (for all walks).
    [free_nb, nr_free_nb] = getFreeNb(walks(:,1:d*c),N,d);
    
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


%% Q.5

% Generate N walks from g(x_(0:n))
N=10000;
n = 10;
d = 2;
walks = zeros(N,d*(n + 1));
weight = 4*ones(N,n);
result = ones(1,n);

for c = 1:n
    
    if c ~= 0
        mult = resampling(weight(:, c));
        matr = multToMatr(mult);
        walks = matr*walks;
    end
    
    % Get which free nb:s there are and how many (for all walks).
    [free_nb, nr_free_nb] = getFreeNb(walks(:,1:d*c),N,d);
    
    % Get new point for each of the N separate walks.
    prev_col = walks(:,d*c - d + 1:d*c);
    new_col = getFreeStep(N,d,prev_col, free_nb, nr_free_nb);
    walks(:,d*c + 1:d*c + d) = new_col;
    
    % calculate weights
    %N_tot(:, c) = nr_free_nb;
    if c>1
        weight(:,c) = nr_free_nb;
    end
    
 
    
end

result = sum(weight,1)/N;

c_cum = ones(1, n);
c_cum(1,1) = result(1,1);
for i = 2:n
    c_cum(1,i) = c_cum(1,i-1)*result(1,i);
end

c_cum
%% Q7
N = 1000;
n = 240;
rep = 100;
cum_sums = zeros(rep, n);

for i = 1:rep
    cum_sums(i, :) = SISR_sample(N, n);
    if mod(i, 5) == 0
        i
    end
end

c_SISR_mean = mean(cum_sums)

nvec = 1:n;
invnvec = 1./nvec;
nroot_cn = c_SISR_mean.^(invnvec);

plot(1:n,nroot_cn)



