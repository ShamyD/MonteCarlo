


% Q 3.

% Generate N walks from g(x_(0:n))
N=1000000;
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
N=1000000;
n=10;
d = 2;
walks = zeros(N,d*(n + 1));
N_tot = 4*ones(N,n);

for c = 1:n
    
    % Get which free nb:s there are and how many (for all walks).
    [free_nb, nr_free_nb] = getFreeNb(walks(:,1:d*c));
    
    % Get new point for each of the N separate walks.
    prev_col = walks(:,d*c - d + 1:d*c);
    new_col = getFreeStep(N,d,free_nb, nr_free_nb);
    walks(:,d*c + 1:d*c + d) = new_col;
    
    %Update nr_free_nb_tot
    N_tot(:, c) = nr_free_nb;
    
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







