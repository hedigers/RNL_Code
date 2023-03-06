
% Percentage Relative Improvement in Average Loss (PRIAL) for different estimators and a specific
% scenario (helper function used in run_sim_tdist)

%% Input:

% p: number of dimensions

% n: number of observations

% mu: the location vector of the multivariate Gaussian (only of interest if lam equal to inf).

% H: true dispersion matrix

% lam: degrees of freedom of the student-t distribution. A value of inf
% corresponds to sampling from a multivariate Gaussian.

% S: number of simulation runs

%% Output:

% Prial_vec: vector of PRIAL values (the size of the vector depends on the number of compared estimators).

%%

function [Prial_vec] = Return_Prial(p, n, mu, H, lam, S)

rng(747, 'twister');

H = SymPDcovmatrix(H);
H_true = H/trace(H)*p;

H_sample_s = NaN(p,p,S);
H_l_s = NaN(p,p,S);
H_nl_s = NaN(p,p,S);
H_rob_l_s = NaN(p,p,S);
H_rob_l_gmv_s = NaN(p,p,S);
H_rob_l_A_s = NaN(p,p,S);
H_rob_l_t_s = NaN(p,p,S);
H_rob_nl_s = NaN(p,p,S);
H_rob_nl_corrbased_s = NaN(p,p,S);

for s = 1:S

    if lam == inf
        X = mvnrnd(mu, H, n);
    else
        sqrts=diag(sqrt(diag(H)).^(-1));
        C=sqrts*H*sqrts;
        X=mvtrnd(C,lam,n);
        X=X.*diag(sqrts).^(-1)';
    end

    meanX=mean(X);
    X=X-meanX(ones(n,1),:);

    H_sample = X'*X/(n-1);
    H_sample = H_sample/trace(H_sample)*p;
    H_sample_s(:,:,s) = H_sample;

    [H_l]=mycov1para(X);
    H_l = H_l/trace(H_l)*p;
    H_l_s(:,:,s) = H_l;

    H_nl = myQIS(X,1);
    H_nl = H_nl/trace(H_nl)*p;
    H_nl_s(:,:,s) = H_nl;

    [H_rob_l] = Robust_Lin_Shrink(X,0);
    H_rob_l = H_rob_l/trace(H_rob_l)*p;
    H_rob_l_s(:,:,s) = H_rob_l;

    [H_rob_l_gmv] = Robust_Lin_Shrink(X,1);
    H_rob_l_gmv = H_rob_l_gmv/trace(H_rob_l_gmv)*p;
    H_rob_l_gmv_s(:,:,s) = H_rob_l_gmv;

    [H_rob_l_A] = Robust_Lin_Shrink(X,2);
    H_rob_l_A = H_rob_l_A/trace(H_rob_l_A)*p;
    H_rob_l_A_s(:,:,s) = H_rob_l_A;

    [H_rob_l_t] = Robust_Lin_Shrink(X,3);
    H_rob_l_t = H_rob_l_t/trace(H_rob_l_t)*p;
    H_rob_l_t_s(:,:,s) = H_rob_l_t;

    H_rob_nl = RNL(X, false, true, true, 1, true);
    H_rob_nl = H_rob_nl/trace(H_rob_nl)*p;
    H_rob_nl_s(:,:,s) = H_rob_nl;

    H_rob_nl_corrbased = RNL(X, true, true, true, 1, true);
    H_rob_nl_corrbased = H_rob_nl_corrbased/trace(H_rob_nl_corrbased)*p;
    H_rob_nl_corrbased_s(:,:,s) = H_rob_nl_corrbased;

end

% Standard PRIAL 
Prial_vec(1) = PRIAL(H_sample_s,H_l_s,H_true);
Prial_vec(2) = PRIAL(H_sample_s,H_nl_s,H_true);
Prial_vec(3) = PRIAL(H_sample_s,H_rob_l_s,H_true);
Prial_vec(4) = PRIAL(H_sample_s,H_rob_l_gmv_s,H_true);
Prial_vec(5) = PRIAL(H_sample_s,H_rob_l_A_s,H_true);
Prial_vec(6) = PRIAL(H_sample_s,H_rob_l_t_s,H_true);
Prial_vec(7) = PRIAL(H_sample_s,H_rob_nl_s,H_true);
Prial_vec(8) = PRIAL(H_sample_s,H_rob_nl_corrbased_s,H_true);

Prial_vec = Prial_vec(:);

end
