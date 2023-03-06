
% Implementation of several robust linear shrinkage estimators

%% Input:

% X: data matrix of dimension nxp, where n is the number of observations
% and p the number of dimensions. p>n is allowed.

% est_type: either 0,1,2 or 3, where 

% 0 (default): the robust linear shrinkage estimator of 
% Y. Chen, A. Wiesel, and A. O. Hero,
% 'Robust shrinkage estimation of high-dimensional covariance matrices',
% IEEE Transactions on Signal Processing, vol. 59, no. 9, pp. 4097–4107, 2011.
%
% 1: the robust linear shrinkage estimator of 
% L. Yang, R. Couillet, and M. R. McKay,
% 'Minimum variance portfolio optimization with robust shrinkage
% covariance estimation', in 2014 48th Asilomar Conference on Signals, Systems and Computers. IEEE, 2014, pp. 1326–1330.
%
% 2: the regularized estimator of 
% T. Zhang and A. Wiesel, 'Automatic diagonal loading for Tyler's
% robust covariance estimator', in 2016 IEEE Statistical Signal Processing Workshop (SSP). IEEE, 2016, pp. 1–5.

% 3: the robust estimator of
% J. Goes, G. Lerman, and B. Nadler, 'Robust sparse covariance estimation by thresholding Tyler's M-estimator',
% Annas of Statistics, vol. 48, no. 1, pp. 86-110, 2020. [Online]. Available: https://doi.org/10.1214/18-AOS1793

%% Output:

% H_hat: a pxp covariance matrix estimate.

%%

function [H_hat] = Robust_Lin_Shrink(X, est_type)

if ~exist('est_type')
    est_type = 0;
end

[n,p] = size(X);

X = X-mean(X);

if est_type == 1

    c_N = p/n;
    eps = 0.1;
    low = eps + max(0,1-(1/c_N));

    [rho_opt, ~] = fminbnd(@(rho) sig_sc(rho,X),low,0.999);
    [H_hat] = C_fixpoint(rho_opt, X, eye(p),1);

elseif est_type == 2

    Y=X./sqrt(sum(X.^2,2));
    Ri = Y'*Y;
    R = 1/n*Ri;

    zeta = p*trace(R^2) - p/n-1;
    rho_opt = 1/n * (zeta + 1 + p)/(zeta + p/n);

    [H_hat] = C_fixpoint(rho_opt, X, eye(p), 2);

elseif est_type == 3

    tresh = sqrt(log(p)/n);
    rho_opt = 10/11;

    [H_alpha] = C_fixpoint(rho_opt, X, eye(p), 1);

    H_tmp = H_alpha - rho_opt*eye(p);
    inner_mat = p*(H_tmp/trace(H_tmp));

    H_hat = zeros(p,p);
    H_hat(abs(inner_mat)>tresh) = inner_mat(abs(inner_mat)>tresh);

else

    X = X./sqrt(sum(X.^2,2));
    Ri = X'*X;
    R = p/n*Ri;

    rho_opt = (p^2 + (1-2/p)*trace(R^2))/((p^2-n*p-2*n) + (n+1 + 2*(n-1)/p)*trace(R^2));

    rho_opt = min((rho_opt>=0)*rho_opt,1);

    [H_hat] = C_fixpoint(rho_opt, X, eye(p), 0);

end

end

