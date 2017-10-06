function valIC = getICvalues_wmm(X, param, labels, cp)

N = size(X,1);

% compte beta
beta_aic = (log(3) - log(log(log(N))))/log(N);
beta_bic = (log(log(N)) - log(log(log(N))))/log(N);
beta_min = (log(log(N)) / (log(N)));
beta_max = 1 - beta_min;
beta_range = [0 beta_aic beta_bic beta_min:(beta_max-beta_min)/20:beta_max];

% compute C_phi_beta(k)
C_phi_beta_k_aic      = (N^beta_aic) * log(log(N));
C_phi_beta_k_bic      = (N^beta_bic) * log(log(N));
C_phi_beta_k_beta_min = (N^beta_min) * log(log(N));
C_phi_beta_k_beta_max = (N^beta_max) * log(log(N));
C_phi_beta_k_beta_range = (N.^beta_range) .* log(log(N));

% Get the information
mu = param.mu;
kappa = param.kappa;
alpha = param.weight;

k = length(alpha);
p = size(mu,2);

% compute P(k)
P_k = (k*(p + 1)) + k - 1;

%
%% Compute log likelihood for the optimal label
c= 0.5;
d = p;
for h=1:k
    K=log(1/chgm(0.5,d/2,kappa(h)));
    %%MAPLE HIGH DIMENSIONAL CALL
    %       [ans]=maple(sprintf('log(1/KummerM(%f,%f,%15.5f))',0.5,d/2,kappa(h)));
    %       K=sscanf(ans,'%f');
    fh_X(:,h)=K+(kappa(h)*(mu(h,:)*X')'.^2);
    ind=fh_X(:,h)>700;
    fh_X(ind,:)=fh_X(ind,:)-1000;
    fh_X(:,h)=exp(fh_X(:,h));
end

llh = sum(log(sum(bsxfun(@times, alpha,fh_X), 2)));


cp = cp + 0.00001;
lcp = log(cp);

sumLCPOp = 0;     % Sum of log conditional probability
for tnumcl=1:k
    tindxs = find(labels==tnumcl);
    sumLCPOp = sumLCPOp + sum(lcp(tindxs, tnumcl));
end

% Compute the IC scores
valIC.BIC = -2*llh + (C_phi_beta_k_bic * P_k);
valIC.AIC = -2*llh + (C_phi_beta_k_aic * P_k);
valIC.ICL = -2*llh + (log(N) * P_k) - (2*sumLCPOp);
valIC.LLHT = llh;
valIC.beta_min = -2*llh + (C_phi_beta_k_beta_min * P_k);
end