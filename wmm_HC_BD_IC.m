% Function STATUS: ACTIVE

% Purpose: Compute Information Criterion score when
% applying hierarchical clustering with Bregman divergence based approach.

% Author: Md. Abul Hasnat

function [allInfo allClust Clusters_Params] = wmm_HC_BD_IC(X, params)

[n,d]=size(X);
a = 0.5;
c = d/2;

% Get the parameters from the structure
k = length(params.weight);
label = params.label;
cp = params.cp;
eta = params.expectation;
theta_cl = params.natural;
weight = params.weight';

% build extension matrices for X and mu : extension matrix
[X_ex M_ex] = getExtensionMatrices(X, params.source.mu, d);

%% Perform hierarchical clustering
numberOfCluster = k;

LNF = zeros(numberOfCluster,1);
DLNF = zeros(numberOfCluster,1);

for j = 1:numberOfCluster
    normTheta(j) = getNormThetaNR(norm(eta(j,:)), d);
    
    % Compute Bragman divergence
    LNF(j) = log(chgm(0.5, d/2, normTheta(j))); % The log normalizing function
    DLNF(j) = (theta_cl(j,:) * eta(j,:)') - LNF(j); % Dual (expectation parameter) of the log normalizing function
end
clear j;

% Compute distance among clusters
indx = 1;
D_G = zeros(numberOfCluster);
Div_G = zeros(1, (numberOfCluster*(numberOfCluster-1))/2);
distMat = zeros(size(Div_G));

for i = 1 : numberOfCluster-1
    for j = i+1 : numberOfCluster
        % Left sided distance
        Div_G(indx) = DLNF(i) - DLNF(j) - ( (eta(i,:) - eta(j,:)) * theta_cl(j,:)');
        distMat(indx) = weight(i) * weight(j) * Div_G(indx);
        indx = indx+1;
    end
end
clear i j;

% Apply hierarchical clustering (use matlab function 'linkage')
Z = linkage(distMat,'average'); % 'average distance'
Clusters_Params = cell(1, numberOfCluster);

%% Perform clustering
for numClust = numberOfCluster : -1 : 1
    ucp = zeros(size(X,1), numClust);
    
    % Do initial clustering
    T = cluster(Z,'maxclust',numClust);
    isConsidered = zeros(size(T));
    
    % Allocate updated cluster parameter vectors
    Up_eta = zeros(numClust, size(eta,2));
    Up_theta_cl = zeros(numClust, size(theta_cl,2));
    Up_weight = zeros(numClust, 1);
    update_mu = zeros(numClust, d);
    update_kappa = zeros(1,numClust);
    
    subset_cl = cell(1, numClust);
    indx=1;

    numClMerge = 0;
    
    % Find the merged clusters and Update cluster centroid (usign Bregman centroid)
    for i=1:numberOfCluster
        if(isConsidered(i)), continue; end

        cl_label = T(i);
        % get the associated subsets of original cluster
        subset_cl{indx} = find(T==cl_label);
        isConsidered(subset_cl{indx}) = 1;
        
        % Merge the subset and update centroid and weight
        if(length(subset_cl{indx}) > 1)
            % Left sided Bregman centroid computation with expectation
            % parameters
            Up_weight(indx) = sum(weight(subset_cl{indx}));
            Up_eta(indx,:) = sum(bsxfun(@times, weight(subset_cl{indx}), eta(subset_cl{indx},:))) ./ Up_weight(indx);
            
            % Convert to source parameters (mu, kappa) from expectation parameter (eta)
            normUp_theta(indx) = getNormThetaNR(norm(Up_eta(indx, :)), d);
            
            % Compute R(normTheta)
            g_norm_theta_Up(indx) = (a/c) * (chgm(a+1, c+1, normUp_theta(indx)) / chgm(a, c, normUp_theta(indx)));
            R_norm_Up_theta(indx) = g_norm_theta_Up(indx) / normUp_theta(indx);
            Up_theta_cl(indx, :) = Up_eta(indx, :) ./ R_norm_Up_theta(indx); % natural parameter
            
            % update source parameters
            update_kappa(indx) = normUp_theta(indx);
            update_W(indx, :) = Up_theta_cl(indx, :) ./ normUp_theta(indx);
            update_mu(indx, :) = getMufromW(update_W(indx, :), d);
            
            % Update conditional probability
            ucp(:, indx) = sum(cp(:, subset_cl{indx}),2);
            
            % Find the number of data samples merged during the merging process
            for ijk = 1:length(subset_cl{indx})
                numClMerge = numClMerge + length(find(label == subset_cl{indx}(ijk)));
            end
        else
            % no update
            Up_weight(indx) = weight(subset_cl{indx});
            Up_eta(indx,:) = eta(subset_cl{indx},:);
            
            % Convert to source parameters (mu, kappa) from expectation parameter (eta)
            normUp_theta(indx) = getNormThetaNR(norm(Up_eta(indx, :)), d);

            % Compute R(normTheta)
            g_norm_theta_Up(indx) = (a/c) * (chgm(a+1, c+1, normUp_theta(indx)) / chgm(a, c, normUp_theta(indx)));
            R_norm_Up_theta(indx) = g_norm_theta_Up(indx) / normUp_theta(indx);
            Up_theta_cl(indx, :) = Up_eta(indx, :) ./ R_norm_Up_theta(indx); % natural parameter

            % update source parameters
            update_kappa(indx) = normUp_theta(indx);
            update_W(indx, :) = Up_theta_cl(indx, :) ./ normUp_theta(indx);
            update_mu(indx, :) = getMufromW(update_W(indx, :), d);
            
            % Update conditional probability
            ucp(:, indx) = cp(:, subset_cl{indx});
        end
        
        indx = indx + 1;
    end
    
    Clusters_Params{numClust}.src_params.mu = update_mu;
    Clusters_Params{numClust}.src_params.kappa = update_kappa;
    Clusters_Params{numClust}.weight = Up_weight;
    Clusters_Params{numClust}.exp_params = Up_eta;
    Clusters_Params{numClust}.nat_params = Up_theta_cl;
    Clusters_Params{numClust}.cp = ucp;
    
    
    %% Compute Information Criterion values
    clear prms
    prms.mu = update_mu;
    prms.kappa = update_kappa;
    prms.weight = Up_weight';
    
    % --> Get the class number for each component
    clear sufStat_minus_expectParam Log_Normalizing_Function gradDLNF innerProdTerm divergence
 
    for j=1:numClust
        normTheta2(j) = getNormThetaNR(norm(Up_eta(j,:)), d);
        
        % Compute R(normTheta)
        g_norm_theta2(j) = (a/c) * (chgm(a+1, c+1, normTheta2(j)) / chgm(a, c, normTheta2(j)));
        R_norm_theta2(j) = g_norm_theta2(j) / normTheta2(j);
        theta_cl2(j, :) = Up_eta(j, :) ./ R_norm_theta2(j); % natural parameter       
        
        LNF2(j) = log(chgm(a, c, normTheta2(j))); % The log normalizing function
        DLNF2(j) = (theta_cl2(j,:) * Up_eta(j,:)') - LNF2(j); % Dual (expectation parameter) of the log normalizing function
        
        GDLNF2(j,:) = theta_cl2(j, :);
        
        sufStat_minus_expectParam(:, :, j) = bsxfun(@minus, X_ex , Up_eta(j, :));
        innerProdTerm(:,j) = sufStat_minus_expectParam(:, :, j) * GDLNF2(j,:)';        
        divergence(:,j) = -(DLNF2(j) + innerProdTerm(:,j));
    end
       
    % Hard clustering
    [~, allClust(:, numClust)] = min(divergence,[], 2);
    
    valIC2 = getICvalues_wmm(X, prms, allClust(:, numClust), ucp);
    
    allInfo.tE(numClust) = -sum(sum(ucp .* log2(ucp)));
    allInfo.numMerge(numClust) = numClMerge;
    
    allInfo.BIC(numClust) = valIC2.BIC;
    allInfo.AIC(numClust) = valIC2.AIC;
    allInfo.ICL(numClust) = valIC2.ICL;
    allInfo.LLH(numClust) = valIC2.LLHT;
    allInfo.Beta_min(numClust) = valIC2.beta_min;
end