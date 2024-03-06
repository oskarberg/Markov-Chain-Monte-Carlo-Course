N = 1000;  % Number of particles
n = 400;     % length of each self-avoiding walk
weights = ones(N, 1);  % Initial weights for all paths

X = zeros(n+1, 2, N);  % Storage for paths (+1 for the origo position)

directions = [0, 1; 0, -1; -1, 0; 1, 0];


cumulativeProduct = 1;

for k = 1:n
    for i = 1:N
        % Get current position
        currentPos = squeeze(X(k,:,i));
        
        % Determine available moves
        availableMoves = [];
        for d = 1:4
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
%cn2_var = var(weights)

% Estimate cn(2) SISR
% How to find the variance?? 
% SISR seems worse than SIS
fprintf('Estimated cn(2) for n=%d: %f\n', n, cumulativeProduct);



