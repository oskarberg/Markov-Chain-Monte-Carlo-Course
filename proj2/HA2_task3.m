clear; close all; 
%%
N = 1000; %Number of particles
n =30; %time
d = 2;
X = zeros(N, d, n); %state space: particles i, state, time k
w = zeros(N, 1, n);  %weights
cn = [];

for i = 1:N   
    w(i, 1, 1) = [1]; 
end


step = [eye(d), -eye(d)]'; 
for k = 2:n+1 
    for i = 1:N  
        number = randi(2*d);   
        X(i, :, k) = step(number, :) + X(i, :, k-1); 
        oneorzero = [];
        for j = 1:k-1    
         indicator = isequal(X(i, :, k), X(i, :, j));    
         oneorzero = [oneorzero, indicator];
        end    
        if ismember(1, oneorzero)    
            w(i, 1, k) = 0;  
        else
            w(i, 1, k) = 4*w(i, 1, k-1);  
        end
    end
    cn = [cn, (1/N)*sum(w(:, :, k))];  
end
mu = cn.^(1./(1:n))
%cn