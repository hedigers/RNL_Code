function [H_hat, rho_opt, sig_opt] = Robust_Lin_Shrink(X, est_type)

% est_type = 0: Robust Shrinkage Estimation of High Dimensional Covariance Estimtors (Chen)
% est_type = 1: GMV (zhang et al., 2014)
% est_type = 2: Regularized M-Estimators of Scatter Matrix (Ollila and Tyler)

% X is nXp
% X is assumed to be centered

if ~exist('est_type')
    est_type=0;
end

[n,p] = size(X);

X=X-mean(X);

if est_type == 1

    %%% Competitor 2: R-GMV-LS: the robust linear shrinkage estimator of (Zhang et al., 2014). This esti-
    %%% mator is especially well suited for global-minimum-variance portfolios.
    %%% There should be no standardization here %%%

    c_N = p/n;
    eps = 0.1;
    low = eps + max(0,1-(1/c_N));

    [rho_opt, sig_opt] = fminbnd(@(rho) sig_sc(rho,X),low,0.999);

    [H_hat] = C_fixpoint(rho_opt, X, eye(p),1);

elseif est_type == 2

    %%% Competitor 3: R-M-LS: the regularized M-estimator of Ollila and Tyler (2014)
    %%% There should be no standardization here %%%

    if n <= p
        beta = (n/p)-(n/(2*p));
        rho = 1 - beta;
    else
        rho = 0;
    end

    [R] = C_fixpoint(rho, X, eye(p), 0);

    rho_opt = (p*trace(R)-1)/(p*trace(R)-1 + n*(p+1)*((1/p)*trace(R^(-2))-1));

    [H_hat] = C_fixpoint(rho_opt, X, eye(p), 1);
    sig_opt = NaN;

elseif est_type == 3
    %%% Competitor 4: R-A-LS: the adapted estimator of Zhang and Wiesel (2016)
    %%% There should be no standardization here %%%

    Y=X./sqrt(sum(X.^2,2));
    Ri = Y'*Y;
    R=1/n*Ri;
    sig_opt = NaN;

    zeta= p*trace(R^2) - p/n-1;
    rho_opt= 1/n * (zeta + 1 + p)/(zeta + p/n);

    [H_hat] = C_fixpoint(rho_opt, X, eye(p), 2);

else

    %%% Competitor 1: R-LS: the robust linear shrinkage estimator of (Chen et al., 2011). An estimator that
    %is "widely used and performs well in practice” (Sun et al., 2014)

    X=X./sqrt(sum(X.^2,2));
    Ri = X'*X;
    R=p/n*Ri;

    rho_opt = (p^2 + (1-2/p)*trace(R^2))/((p^2-n*p-2*n) + (n+1 + 2*(n-1)/p)*trace(R^2));

    rho_opt = min((rho_opt>=0)*rho_opt,1);
    sig_opt = NaN;

    [H_hat] = C_fixpoint(rho_opt, X, eye(p), 0);

end

end

