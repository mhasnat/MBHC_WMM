% Generate samples from Watson distribution
% Used the MoVMF code was written by Arindam Banerjee and Suvrit Sra.
% Author: Md. Abul Hasnat
% Date: 10/6/2013
function [wmmSample, labels, mixture] = Samples_Watson_mixture_model_gui(num_comp, num_samp)
b_s_eq_s = 0;

dim = 3;

% Set the dimensionality ... corr to 'd'
mixture.dim = dim;

% Equal priors for all clusters.
% if both sides are not equally spaced then
if(~b_s_eq_s)
    mixture.priors = ones(1,(num_comp*2))/(num_comp);
    p = rand(1,num_comp);
    pm1 = abs(1-p);
    mixture.priors(1:2:end) = mixture.priors(1:2:end) .* p;
    mixture.priors(2:2:end) = mixture.priors(2:2:end) .* pm1;
else
    mixture.priors = ones(1,(num_comp*2))/(num_comp*2);
end

% Record no. of clusters we are dealing with.
mixture.num_clus = num_comp*2;

% Set the centroids
mixture.centers = zeros(num_comp*2,dim);

mixture.centers(1,:) = unitrand(dim); % Now setting it randomly
mixture.centers(2,:) = -mixture.centers(1,:); % Vector in opposite dimension

for i = 2:num_comp   
    while(1)
        sampleRand = unitrand(dim); % Now setting it randomly
        % The sample center should be far enough (>1/num_comp) from the other samples
        distAmongCenters = 1-abs(mixture.centers(1:(i*2),:) * sampleRand);
        if(min(distAmongCenters)>1/(2*num_comp-1))
            mixture.centers((i*2)-1,:) = unitrand(dim); % Now setting it randomly
			mixture.centers((i*2),:) = -mixture.centers((i*2)-1,:); % Vector in opposite dimension
            break;
        end
    end
end

% The covars(:,:,i) corresponds to \Sigma_i for the mixture.
tkappas = 30 + randi(50,num_comp,1); % Now setting it randomly
mixture.kappas(1:2:num_comp*2) = tkappas;
mixture.kappas(2:2:num_comp*2) = tkappas;

%% Generate samples
% mixture.centers
% mixture.kappas'
[wmmSample,labels] = emsamp(mixture, num_samp);