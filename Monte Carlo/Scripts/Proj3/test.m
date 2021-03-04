load('C:\Users\Seamu\MATLAB\Projects\MonteCarlo\Monte Carlo\Data\coal_mine_disasters.mat')
plot(T)

%%

N = 10^4*1.25; % count in burn in time

d = 3;%amount of breakpoints
t1 = 1658;
td1 = 1980;

t=zeros(N, d+1);
t(:,1) = t1;
t(:,end) = td1;

%%
    a = zeros(10^3, 4);
    for i=1:10^3
        a(i,:) = rwp([1 3 4 7], 0.5);
    end
    mean(a)