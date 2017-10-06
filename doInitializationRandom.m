function [clust mu] = doInitializationRandom(X, C)
[n,d]=size(X);
%Initialize Means
for i=1:C
    mu(i,:) =-1+2*rand(1,d);
    %mu(i,:) = mean(X)+ 0.01*rand(1,d);
    mu(i,:)=mu(i,:)/norm(mu(i,:)+1e-10);
end

% assign points to nearest cluster
simMat        =  (X*mu').^2;
[simax,clust] =  max(simMat,[],2);
