
% Updating the eigenvectors

%% Input:

% Z: data matrix of dimension nxp, where n is the number of observations
% and p the number of dimensions. p>n is allowed.

% eig_fix: a p-dimensional vector of eigenvalues (returned from QIS).

% alpha: value between 0 and 1 (default is 1), representing the convex
% combination between the sample eigenvalues and the fixed eigenvalues of
% QIS (within VIteration).

% crit: either 0 or 1 (default is 1), specifying the convergence criterion.

% H_start: initial covariance matrix estimate. Default is the identity matrix.

%% Output:

% H_hat: a pxp covariance matrix estimate.

% V: a pxp matrix of eigenvectors.

%%

function [H_hat, V] = VIteration(Z, eig_fix, alpha, crit, H_start)

tol = 1e-5;
maxit = 10;
[n,p] = size(Z);
V = eye(p);

eig_fix_inv = eig_fix.^(-1);

i = 0;
diagnorm = 1000;

if nargin < 5
    lam = ones([p,1]);
    H_new = V*((1-alpha)*diag(lam) + alpha*diag(eig_fix))*V';
else
    H_new = H_start;
end

lower = diag((Z*inv(H_new)*Z')./p);
ratioV = ((Z./lower)'*(Z))./n;
ratioV = (ratioV+ratioV')/2;

while abs(diagnorm) > tol && i < maxit

    i = i+1;

    V_old = V;
    ratioV_old = ratioV;
    
    H_old = H_new;
    Hinvold = inv(H_old);

    lower = diag((Z*Hinvold*Z')./p);
    ratioV = ((Z./lower)'*(Z))./n;
    ratioV = (ratioV+ratioV')/2;

    [V,lam] = eig(ratioV,'vector');
    [~, ind] = sort(lam);
    V = V(:,ind);
    lam = lam(ind);

    H_new = V*((1-alpha)*diag(lam) + alpha*diag(eig_fix))*V';

    if crit == 1 && alpha ~= 1
        diagnorm = norm(H_new-H_old);
    else
        diagnorm = norm(V_old'*ratioV_old*V_old*diag(eig_fix_inv) - diag(eig_fix_inv)*V'*ratioV*V, 'fro');
    end

end

H_hat = H_new;

end