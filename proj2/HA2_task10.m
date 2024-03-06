clear all
load population_2024.mat
%% 
for i = 1:1000
N = 10000;
n = 100;
A = 0.8;
B = 3.8;
C = 0.6;
D = 0.99;
G = 0.8; 
H = 1.25;
tau = zeros(1,n+1); % vector of filter means
w = zeros(N,1);
p = @(x, y) (y >= G*x).*(y <= H*x)./((H-G).*x); % observation density, for weights
part = unifrnd(C, D, N, 1); % initialization
w = p(part, Y(1)); % weighting
tau(1) = sum(part.*w)/sum(w); % estimation

% conf int
[xx,I]=sort(part); % sort data
cw=cumsum(w(I))/sum(w); % cumulative normalized weightsum
% for sorted data
Ilower=find(cw>=0.025,1); % index for lower 2.5% quantile
Iupper=find(cw>=0.975,1); % index upper 2.5% quantile
taulower(1)=xx(Ilower); % lower 2.5% quantile
tauupper(1)=xx(Iupper); % upper 2.5% quantile

ind = randsample(N,N,true,w); % selection
part = part(ind);

for k = 1:n % main loop 
    R = unifrnd(A, B, N, 1);
    part = R.*part.*(1-part); % mutation
    w = p(part, Y(k+1)); % weighting
    tau(k + 1) = sum(part.*w)/sum(w); % estimation
    
    % conf int
    [xx,I]=sort(part); % sort data
    cw=cumsum(w(I))/sum(w); % cumulative normalized weightsum
    % for sorted data
    Ilower=find(cw>=0.025,1); % index for lower 2.5% quantile
    Iupper=find(cw>=0.975,1); % index upper 2.5% quantile
    taulower(k+1)=xx(Ilower); % lower 2.5% quantile
    tauupper(k+1)=xx(Iupper); % upper 2.5% quantile
     
    ind = randsample(N,N,true,w); % selection
    part = part(ind);
end

% plot(tau)
% hold on
% plot(X)
% hold off
% legend('Tau', 'Real X');
% title('Tau vs Real X N = 10000')

count = 0; %start counter
index = 1; %index of which X(j) is out of the interval
out = []; 
for j = 1:length(X)
    if X(j) <= tauupper(j) && X(j) >= taulower(j)
        count = count+1;
    else
        out(index) = j;
        index = index+1;
    end
end

x_tau(i) = mean(X-tau');
var_x(i) = var(X-tau');
ci_width(i) = mean(tauupper-taulower);
out_ci(i) = length(out);
end
mean(x_tau)
mean(var_x)
mean(ci_width)
mean(out_ci)