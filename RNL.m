
% Robust nonlinear shrinkage

%% Input:

% X: data matrix of dimension nxp, where n is the number of observations
% and p the number of dimensions. p>n is allowed.

% complicated: either 0 or 1 (default is 1). If equal to 0, VIteration is not called (the
% eigenvectors of the sample covariance matrix are taken).

% normscale: either 0 or 1 (default is 1). If equal to 1, the data is scaled by the Euclidean norm of X.  

% corrbased: either 0 or 1 (default is 1). If equal to 1, the data is scaled by the p sample standard deviations of X. 

% beta: value between 0 and 1 (default is 1), representing the convex
% combination between the sample eigenvalues and the fixed eigenvalues of
% QIS (within VIteration).

% addQIS: either 0 or 1 (default is 1). If equal to 1, the eigenvalues are
% updated in the end using the optimized eigenvectors of VIteration. If equal to 0, QIS is only once applied in the
% beginning to get the fixed eigenvalues.

%% Output:

% Sigma: a pxp covariance matrix estimator.

%% Examples:

% R-C-NL corresponds to calling RNL(X)

% R-NL corresponds to calling RNL(X,0)

% SRTy of A. Breloy, E. Ollila, and F. Pascal, 'Spectral shrinkage of Tyler's M-estimator of covariance
% matrix' corresponds to calling RNL(X,1,1,0,beta,0), where the user has
% to choose a value for beta = alpha/(1+alpha).

%%

function [Sigma] = RNL(X, corrbased, complicated, normscale, beta, addQIS)

if nargin < 2
    complicated = 1;
    normscale = 1;
    corrbased = 1;
    beta = 1;
    addQIS = 1;
end

if nargin == 2
    complicated = 1;
    normscale = 1;
    beta = 1;
    addQIS = 1;
end

[n,p] = size(X);

location = mean(X);
x = X-location;

dsquare = var(x);
if corrbased == 0
    dsquare = ones(1,p);
end

x = x./sqrt(dsquare);
if normscale == 1
    Z = x./sqrt(sum(x.^2,2));
    [~,eig_fix] = myQIS(Z,1);
else
    [~,eig_fix] = myQIS(x,1);
    Z = x;
end

[eig_fix] = sortLambda(eig_fix, n, p);

if complicated == 1
    [H_tmp,V] = VIteration(Z,eig_fix,beta,1);
else
    [V,lam] = eig(cov(Z),'vector');
    [~, ind] = sort(lam);
    V = V(:,ind);
    H_tmp = V*diag(eig_fix)*V';
end

if beta == 1 && addQIS == 1
    Hinv = V*diag(eig_fix.^(-1))*V';
    lower = diag((Z*Hinv*Z')./p);
    Y = (Z./sqrt(lower));
    H0 = myQIS(Y,1);
else
    H0 = H_tmp;
end

Sigma = diag(sqrt(dsquare))*H0*diag(sqrt(dsquare));
Sigma = SymPDcovmatrix(Sigma);

Sigma = trace(cov(X)).*Sigma./trace(Sigma);

end

