function [ClustFinal] = movmf_with_Initialization(vectors,k, initClust)

[D,V] = size(vectors);
dim   = V;

mu = initClust.mu;
clust = initClust.clust;

Clust1 = clust;

%----------------------------------------------
% You can cut the code at this point
%----------------------------------------------

% --------------------------------
% movMF iterations
% --------------------------------

diff      = 1;
epsilon   = 0.0001;
value     = 100;
iteration = 1;

kappaMax = 5000;
kappaMin = 1;

% initializing kappa, alpha
kappa    = kappaMin*ones(1,k);
for h = 1:k
  alpha(h) = length(find(clust==h))/D;
end

T = 1;
t = 1;

% display('Starting main iterations ...');
while (diff > epsilon)
%while (iteration < 10)

%   display(['Iteration ',num2str(iteration)]);
  iteration = iteration + 1;
  oldvalue   = value;
  
  % assignment of points to nearest vMFs
  
  logNormalize  = log(alpha) + (dim/2-1)*log(kappa) - (dim/2)*log(2*pi) - logbesseli(dim/2-1,kappa); 
  logProbMat    = (vectors*(mu'.*(ones(dim,1)*kappa)) + ones(D,1)* ...
      logNormalize)/T; 
  lpm = logProbMat;
  logSum        = log(sum(exp(logProbMat),2)); % this step, without
                                               % the 1/T in the
                                               % previous line, leads to
                                               % Inf. We do it w/o
                                               % 1/T in C++ using NTL
					       
  logProbMat    = logProbMat - logSum*ones(1,k);
  
  value = sum(sum(logProbMat));
  
  
  % updating component parameters
  alpha  = sum(exp(logProbMat)); % this step leads to Inf
  mu     = exp(logProbMat')*vectors;
  
  for h=1:k
    oldkappa(h) = kappa(h);

    normMu   = sqrt(mu(h,:)*mu(h,:)');
    rbar(h)  = normMu/(t*alpha(h));
    
    mu(h,:)  = mu(h,:)/normMu;
    kappa(h) = (rbar(h)*dim - rbar(h)^3)/(1-rbar(h)^2);
    alpha(h) = alpha(h)/D;
%     display(['Size = ',num2str(D*alpha(h)),' rbar =',num2str(rbar(h)),' kappa = ',num2str(kappa(h))]);
  end

  diff = abs(value - oldvalue);
  
end


[simax,clust] = max(logProbMat,[],2);  
% subplot(2,1,2),plot(1:D,clust,'ro');
ClustFinal = clust;