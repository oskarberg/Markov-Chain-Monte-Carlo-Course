N = 1000;  % Number of particles
n = 50;     % length of each self-avoiding walk

%cn2 = zeros(200,1)
%for n = 1: 200

weights = ones(N, 1);  % Initial weights for all paths

X = zeros(n+1, 2, N);  % Storage for paths (+1 for origo position)

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
            %weights(i) = 1*size(availableMoves,1);            %SISR
            weights(i) = weights(i)*size(availableMoves,1);  %SIS
        end
    end
    
   
    % %Resampling step
    % normalizedWeights = weights / sum(weights);  
    % indices = randsample(N, N, true, normalizedWeights);
    % X = X(:,:,indices);
    % cumulativeProduct = cumulativeProduct * mean(weights);


end


% cn2(n) = mean(weights)^(1/n);
% cn2(n)
% 
% end

cn2 = mean(weights)^(1/n)


% Estimate cn(2) SIS
% cn2 = mean(weights)
% cn2_var = var(weights)
% fprintf('Estimated cn(2) for n=%d: %f\n', n, cn2);