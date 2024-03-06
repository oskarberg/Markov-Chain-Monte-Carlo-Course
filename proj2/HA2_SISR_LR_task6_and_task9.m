N = 1000;  % Number of particles
%n = 10;     % length of each self-avoiding walk

dim = 2;
directions = [eye(dim); -eye(dim)];

save_beta_est= zeros(3,20);
for run=1:20

cn = zeros(1,15);

for n=1:15

cumulativeProduct = 1;
weights = ones(N, 1);  % Initial weights for all paths

X = zeros(n+1, dim, N);  % Storage for paths (+1 for the origo position)

for k = 1:n
    for i = 1:N
        % Get current position
        currentPos = squeeze(X(k,:,i));
        
        % Determine available moves
        availableMoves = [];
        for d = 1:size(directions,1)
            nextPos = currentPos + directions(d,:);
            if ~any(ismember(squeeze(X(1:k,:,i)), nextPos, 'rows'))
                availableMoves = [availableMoves; directions(d,:)];
            end
        end
        
        % Check if dead end
        if isempty(availableMoves)
            weights(i) = 0; % Particle has reached a dead end
        else
            % Select a move randomly
            moveIdx = randi(size(availableMoves, 1));
            X(k+1,:,i) = currentPos + availableMoves(moveIdx,:);
            weights(i) = 1*size(availableMoves,1);            %SISR
            %weights(i) = weights(i)*size(availableMoves,1);  %SIS
        end
    end
    
   
    %Resampling step
    normalizedWeights = weights / sum(weights);  
    indices = randsample(N, N, true, normalizedWeights);
    X = X(:,:,indices);
    cumulativeProduct = cumulativeProduct * mean(weights);

    % [xx, I ] = sort(weights);
    % cw=cumsum(weights(I))/sum(weights);
    % 
    % Ilower=find(cw>=0.025,1); % index for lower 2.5% quantile
    % Iupper=find(cw>=0.975,1); % index upper 2.5% quantile
    % taulower(k+1)=xx(Ilower); % lower 2.5% quantile
    % tauupper(k+1)=xx(Iupper); % upper 2.5% quantile

end


% Estimate cn(2) SIS
%cn2 = mean(weights);
%var(weights)

% Estimate cn(2) SISR
% How to find the variance?? seems higher than SIS
fprintf('Estimated cn(2) for n=%d: %f\n', n, cumulativeProduct);

cn(n) = cumulativeProduct;%^(1/n)
%cn(n)

end


% Combine into a single matrix X
X = [ones(n, 1) (1:n)' log(1:n)'];

Y = log(cn)'; %+ log((1:n)')

beta_est = X\Y;

beta_est_transformed = zeros(size(beta_est));
beta_est_transformed(1:2) = exp(beta_est(1:2));
beta_est_transformed(3) = beta_est(3) +1;
save_beta_est(:,run) = beta_est_transformed


% Use linear regression variance estimates?


end

Y_pred = X * beta_est; % Predicted values
residuals = Y - Y_pred; % Errors (Residuals)

var_errors = var(residuals)

% Residual diagnostics. Gaussian distributed.
qqplot(residuals)

beta_cov = X'*X*var_errors


var(save_beta_est')'
mean(save_beta_est')'

% 
% fitted = X*beta_est;
% 
% plot(fitted)
% hold on
% plot(Y)


%% Standardized residuals vs fitted

% Calculate residuals
%residuals = Y - X*beta_est;
fitted = X*beta_est;

% Calculate the standard deviation of the residuals
sigma = std(residuals);

% Standardize residuals
std_residuals = residuals / sigma;


% Calculate the percentage of standardized residuals within +/-1.96
%withinBounds = sum(abs(std_residuals) <= 1.96) / length(std_residuals) * 100;

% Plot standardized residuals
figure;
scatter(fitted, std_residuals);
hold on;
line([1, n+1], [1.96, 1.96], 'Color', 'red', 'LineStyle', '--');
line([1, n+1], [-1.96, -1.96], 'Color', 'red', 'LineStyle', '--');
xlabel('Fitted');
ylabel('Standardized Residuals');


%% Graham's asymptotic bound
dim = [2, 3, 5, 10];
Graham_mu =  2*dim - 1 - 1./(2*dim).^2 - 16./(2*dim).^3 
