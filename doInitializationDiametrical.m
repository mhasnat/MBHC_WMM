function [clust mu] = doInitializationDiametrical(X, C)
[n,d]=size(X);
%INITIALIZATION
% --------------------------------
% Getting the initial random means
% --------------------------------
% getting global mean
sumv  = sum(X);
normv = sqrt(sumv*sumv');
mu0   = sumv./normv;

% perturbing global mean to get initial cluster centroids
perturb = 0.1;
for h = 1:C
    randVec     = rand(1,d) - 0.5;
    randNorm    = perturb*rand;
    
    smallRandVec =  randNorm*randVec/sqrt(randVec*randVec');
    mu(h,:)      =  mu0 + smallRandVec;
    mu(h,:)      =  mu(h,:)/sqrt(mu(h,:)*mu(h,:)');
end

% --------------------------------
% Getting the means from spkmeans
% --------------------------------
diff    = 1;
epsilon = 0.001;
value   = 100;
iteration = 1;
while (diff > epsilon && iteration<200)
    
    %   display(['Iteration ',num2str(iteration)]);
    iteration = iteration+1;
    oldvalue      = value;
    
    % assign points to nearest cluster
    simMat        =  (X*mu').^2;
    [simax,clust] =  max(simMat,[],2);
    
    % compute objective function value
    value         = sum(simax);
    
    % compute cluster centroids
    clust_cnt = 1;
    for h=1:C
        clustData = X(find(clust==h),:);
        [~,~,V] = svds(clustData,1);
        if(isempty(V)) continue; end
        mu(clust_cnt,:) = V';
        clust_cnt = clust_cnt + 1;
    end
end