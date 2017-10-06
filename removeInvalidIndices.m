function params = removeInvalidIndices(params)

indx = [];

for i=1:length(params.weight)
    % check alpha
    if(params.weight(i)==0)
        indx = [indx i];
        continue;
    end
    
    % check mu
    mu = params.source.mu(i,:);
    
    if(isnan(mu))
        indx = [indx i];
        continue;
    end
end

% Remove invalid param
params.expectation(indx, :) = [];
params.natural(indx, :) = [];

params.source.kappa(indx) = [];
params.source.mu(indx, :) = [];

params.weight(indx) = [];
params.cp(:, indx) = [];
