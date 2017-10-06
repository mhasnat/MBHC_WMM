function [Init params Rparams logLikelihood] = getWMMParams_bd(X, C, mtype)

% X   -- data (rows=data,cols=dimensions)
% C       -- number of classes
% kappa   -- parameter (distribution)


[n,d]=size(X);

clust = []; mu = [];

if(strcmp(mtype, 'dm'))
    [clust mu] = doInitializationDiametrical(X, C);
elseif(strcmp(mtype, 'rnd'))
    [clust mu] = doInitializationRandom(X, C);
end
% save 'clustMu.mat' 'clust' 'mu'
% load 'clustMu.mat'
Init.label = clust;
Init.mu = mu;

% build extension matrices for X and mu : extension matrix
[X_ex M_ex] = getExtensionMatrices(X, mu, d);

% initializing kappa, alpha
kappaMin = 10;
kappa    = kappaMin*ones(1,C);
a = 0.5;
c = d/2;

for j = 1:C
  indices = find(clust==j);  
  alpha(j) = length(indices)/n;
  
  % parameters W (from mu) and D (from x)
  eta(j,:) = mean(X_ex(indices, :)); 
  normTheta(j) = getNormThetaNR(norm(eta(j,:)), d);
  
  % Compute R(normTheta)
  g_norm_theta(j) = (a/c) * (chgm(a+1, c+1, normTheta(j)) / chgm(a, c, normTheta(j)));
  R_norm_theta(j) = g_norm_theta(j) / normTheta(j);
  theta_cl(j, :) = eta(j, :) ./ R_norm_theta(j); % natural parameter
  
  % compute V = W'*eta
  %   W(j,:) = M_ex(j, :);
  %   V(j) = W(j,:) * eta(j,:)';
  
  % Compute kappa using Newton-Raphson method
  %   kappa(j) = getKappausingNR(V(j), d);
  
  % Compute theta
  %   theta_cl(j,:) = kappa(j) * W(j,:);
  
  % % %   % The log normalizing function
  % % %   LNF(j) = log(chgm(0.5, d/2, kappa(j)));
  % % %
  % % %   % Dual (expectation parameter) of the log normalizing function
  % % %   DLNF(j) = (theta_cl(j,:) * eta(j,:)') - LNF(j);
end
oldParams.alpha = alpha;
oldParams.eta = eta;
oldParams.kappa = normTheta;
oldParams.mu = mu;
oldParams.theta_cl = theta_cl;
% showWatsonClustering(X, clust, size(mu,1), mu);

fh_X=zeros(n,C);
ph_X=zeros(n,C);
K=0;

%EM ALGORITHM
cnt=1;
MAX_ITERATIONS = 200;

kappa_old=kappa;
mu_old=mu;

llh_old = 0;
llh_new = 100;

diffLLH = abs(llh_old - llh_new);
llh_Th = 0.001;

% sorce parameters
update_kappa = normTheta;
update_mu = mu;

while(cnt<MAX_ITERATIONS && diffLLH > llh_Th)
    cnt=cnt+1;

    %EXPECTATION STEP
    llh_old = llh_new;
    
    % Expectation Step
    for j=1:C
        % Compute Bragman divergence
        LNF(j) = log(chgm(0.5, d/2, normTheta(j))); % The log normalizing function
        DLNF(j) = (theta_cl(j,:) * eta(j,:)') - LNF(j); % Dual (expectation parameter) of the log normalizing function
        
        GDLNF(j,:) = theta_cl(j, :);
        
        sufStat_minus_expectParam(:, :, j) = bsxfun(@minus, X_ex , eta(j, :));
        innerProdTerm(:,j) = sufStat_minus_expectParam(:, :, j) * GDLNF(j,:)';
        
        % Compute divergence
        divergence(:,j) = -(DLNF(j) + innerProdTerm(:,j));
        expTerm(:,j) = exp(-divergence(:,j));
        probTerm(:, j) = alpha(j) * expTerm(:,j);
    end

    probTerm = bsxfun(@rdivide, probTerm, sum(probTerm, 2));
    [~, clust] = max(probTerm,[], 2);
    
    % Maximization Step
    sumProbTerm = sum(probTerm);
    
    % Update weight
    alpha = sumProbTerm ./ size(probTerm,1);
    
    % Update parameter
    for j=1:C
        eta(j, :) = sum(bsxfun(@times, probTerm(:,j) , X_ex)) ./ sumProbTerm(j);
        
        % Convert to source parameters (mu, kappa) from expectation parameter (eta)
        normTheta(j) = getNormThetaNR(norm(eta(j,:)), d);
  
        % Compute R(normTheta)
        g_norm_theta(j) = (a/c) * (chgm(a+1, c+1, normTheta(j)) / chgm(a, c, normTheta(j)));
        R_norm_theta(j) = g_norm_theta(j) / normTheta(j);
        theta_cl(j, :) = eta(j, :) ./ R_norm_theta(j); % natural parameter
        
        % update sorce parameters
        update_kappa(j) = normTheta(j);
        update_W(j, :) = theta_cl(j, :) ./ normTheta(j);
        update_mu(j, :) = getMufromW(update_W(j, :), d);
    end
    
    % Compute log likelihood
    for j=1:C
        K=log(1/chgm(0.5,d/2,update_kappa(j)));
        fh_X(:,j)=K+(update_kappa(j)*(update_mu(j,:)*X')'.^2);
        ind=fh_X(:,j)>700;
        fh_X(ind,:)=fh_X(ind,:)-1000;
        fh_X(:,j)=exp(fh_X(:,j));
    end
    logLikelihood = sum(log(sum(bsxfun(@times, alpha,fh_X), 2)));

    llh_new = logLikelihood;
    diffLLH = abs(llh_old - llh_new);
    
    kappa_old=update_kappa;
    mu_old=update_mu;
end

cp = probTerm;

%% Save Parameters
params.expectation = eta;
params.natural = theta_cl;
params.source.kappa = update_kappa;
params.source.mu = update_mu;
params.weight = alpha;
params.cp = probTerm;
params.label = clust;

% Refined parameters by HARD clustering
% Why? Many sample has little probability to become a member of any cluster. 
% Therefore when we add these conditional probabilities for any certain
% cluster, then it shows that it has prior probabiliy > 0. On the other
% hand when we do hard clustering and compute the class conditional
% probability (CP) after that, then we observe that CP == 0, means there is
% no sample that belongs to that cluster.

[~, clustH] = max(cp,[], 2);
numSamp = length(clustH);
cnt = 1;
empIndx = [];
for i=1:length(params.weight)
    prbCl = length(find(clustH==i))/ numSamp;
    if(prbCl>0)
        weight(cnt) = prbCl;
        cnt = cnt + 1;
    else
        empIndx = [empIndx i];
    end
end

Rparams = params;
Rparams.expectation(empIndx,:) = [];
Rparams.natural(empIndx,:) = [];
Rparams.source.kappa(empIndx) = [];
Rparams.source.mu(empIndx,:) = [];
Rparams.weight = weight;

Rcp = cp;
Rcp(:,empIndx) = [];
[~, clustR] = max(Rcp,[], 2);

Rparams.cp = Rcp;
Rparams.label = clustR;