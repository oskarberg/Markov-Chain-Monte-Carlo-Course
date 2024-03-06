%% HA3

close all
clear 
load('coal_mine_disasters.mat')

num_breakpoints = 5;
d = num_breakpoints + 1;
t = linspace(1658, 1980, d+1);

switch num_breakpoints
    case 1
        rho = 0.01;
    case 2
        rho = [0.005 0.017];
    case 3
        rho = [0.024 0.006 0.009];
    case 4
        rho = [0.075 0.055 0.01 0.010];
    case 5
        %rho = 2*[1 1 1 1 1];
        %rho = 0.001*[1 1 1 1 1];
        rho = [0.075 0.055 0.01 0.053 0.025];
    case 6
        rho = [0.08 0.08 0.07 0.020 0.020 0.025];
    otherwise
        error('Unsupported number of breakpoints.');
end


%rho = [0.005 0.017];
%rho = [0.024 0.006 0.009];
%rho = [0.075 0.055 0.01 0.010];
%rho = [0.075 0.055 0.01 0.053 0.025];
%rho = [0.08 0.07 0.05 0.020 0.020 0.025];

psi = 3;
theta = gamrnd(2, 1/psi);
lambda = gamrnd(2, 1/theta, 1, d);

burn_in = 1e3;
for i = 1: burn_in

    theta = drawTheta(lambda, psi);
    lambda = drawLambda(theta, t, tau);
    [~,t] = drawBP(lambda, t, tau, rho);

end

samples = 1e5;

t_save = zeros(samples, length(t));
theta_save = zeros(samples, 1);
lambda_save = zeros(samples, d);
accepted_save = false(samples, num_breakpoints);

for i = 1:samples

    theta = drawTheta(lambda, psi);
    lambda = drawLambda(theta, t, tau);
    [accepted,t] = drawBP(lambda, t, tau, rho);
    
    t_save(i, :) = t;
    theta_save(i) =  theta;
    lambda_save(i,:) = lambda;
    accepted_save(i,:) = accepted; 
end

figure(1)
plot(t_save(:, 2:end-1))
title('Trace Plots')
xlabel('Sample')
ylabel('Year')
ylim([1658 1980])



Acceptance_rate = mean(accepted_save)


figure(2)
% Plot the actual disaster data
plot(tau, 1:length(tau), 'o');
hold on;

startPoint = 0;
cumulativeIntensity = 0;

mean_lambda = mean(lambda_save);
for i = 1:length(mean_lambda)
    % Calculate the endpoint of the current segment
    endPoint = cumulativeIntensity + mean_lambda(i) * (t(i+1) - t(i));
    
    % Plot the line segment
    line([t(i) t(i+1)], [startPoint endPoint], 'Color', [rand rand rand], 'LineWidth', 2);
    
    % Update the starting point for the next segment
    startPoint = endPoint;
    cumulativeIntensity = endPoint;
end


for i = 2:length(t)-1
    line([t(i) t(i)], [0 startPoint], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
end

title('Total number of accidents and expected intensity during 1658-1980');
xlabel('Year');
ylabel('Cumulative number of accidents/expected intensity');

%axis([tau(1) tau(end) 0 max(1:length(tau), startPoint)]);
%legend('Actual Data', 'Expected Intensity', 'Breakpoints', 'Location', 'Best');
hold off;



%% Histogram

figure(2)
for i = 2:d

    histogram(t_save(:,i), 20)
    hold on
end
xlim([1710 1980])

%% Covariance

% Assuming 'samples' contains your MCMC samples
lags = 100; % Number of lags

for i = 2:d
    figure
    autocorr(t_save(:, i), lags);
end

