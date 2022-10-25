function [Prial_vec] = Return_Prial(p, n, mu, H, lam, S)

H = SymPDcovmatrix(H);
H_true = H/trace(H)*p;

H_sample_s = NaN(p,p,S);
H_l_s = NaN(p,p,S);
H_nl_s = NaN(p,p,S);
H_rob_l_s = NaN(p,p,S);
H_rob_l_gmv_s = NaN(p,p,S);
H_rob_l_A_s = NaN(p,p,S);
H_rob_nl_s = NaN(p,p,S);
H_rob_nl_corrbased_s = NaN(p,p,S);

for s = 1:S

    if lam == -1
        X = mvnrnd(mu, H, n);
    else
        sqrts=diag(sqrt(diag(H)).^(-1));
        C=sqrts*H*sqrts;
        X=mvtrnd(C,lam,n);
        X=X.*diag(sqrts).^(-1)';
    end

    meanX=mean(X);
    X=X-meanX(ones(n,1),:);

    H_sample = X'*X/(n-1); % cov(X);
    H_sample = H_sample/trace(H_sample)*p;
    H_sample_s(:,:,s) = H_sample;

    [H_l,shrinkage]=cov1para(X);
    H_l = H_l/trace(H_l)*p;
    H_l_s(:,:,s) = H_l;

    H_nl = QIS(X,1);
    H_nl = H_nl/trace(H_nl)*p;
    H_nl_s(:,:,s) = H_nl;

    [H_rob_l, rho_opt, sig_opt] = Robust_Lin_Shrink(X,0);
    H_rob_l = H_rob_l/trace(H_rob_l)*p;
    H_rob_l_s(:,:,s) = H_rob_l;

    [H_rob_l_gmv, rho_opt, sig_opt] = Robust_Lin_Shrink(X,1);
    H_rob_l_gmv = H_rob_l_gmv/trace(H_rob_l_gmv)*p;
    H_rob_l_gmv_s(:,:,s) = H_rob_l_gmv;

    [H_rob_l_A, rho_opt, sig_opt] = Robust_Lin_Shrink(X,3);
    H_rob_l_A = H_rob_l_A/trace(H_rob_l_A)*p;
    H_rob_l_A_s(:,:,s) = H_rob_l_A;

    H_rob_nl = RNL(X, true, true, false);
    H_rob_nl = H_rob_nl/trace(H_rob_nl)*p;
    H_rob_nl_s(:,:,s) = H_rob_nl;

    H_rob_nl_corrbased = RNL(X, true, true, true);
    H_rob_nl_corrbased = H_rob_nl_corrbased/trace(H_rob_nl_corrbased)*p;
    H_rob_nl_corrbased_s(:,:,s) = H_rob_nl_corrbased;

end

% Standard PRIAL 
Prial_vec(1) = PRIAL(H_sample_s,H_l_s,H_true);
Prial_vec(2) = PRIAL(H_sample_s,H_nl_s,H_true);
Prial_vec(3) = PRIAL(H_sample_s,H_rob_l_s,H_true);
Prial_vec(4) = PRIAL(H_sample_s,H_rob_l_gmv_s,H_true);
Prial_vec(5) = PRIAL(H_sample_s,H_rob_l_A_s,H_true);
Prial_vec(6) = PRIAL(H_sample_s,H_rob_nl_s,H_true);
Prial_vec(7) = PRIAL(H_sample_s,H_rob_nl_corrbased_s,H_true);

Prial_vec = Prial_vec(:);

% disp(['Linear Shrinkage: ', num2str(Prial_vec(1))])
% disp(['QIS: ', num2str(Prial_vec(2))])
% disp(['QuEST: ', num2str(Prial_vec(3))])
% disp(['Robust Linear Shrinkage: ', num2str(Prial_vec(4))])
% disp(['Robust GMV Linear Shrinkage: ', num2str(Prial_vec(5))])
% disp(['Robust M Linear Shrinkage: ', num2str(Prial_vec(6))])
% disp(['RQIS: ', num2str(Prial_vec(7))])
% disp(['RQuEST: ', num2str(Prial_vec(8))])
% disp(['RCQIS: ', num2str(Prial_vec(9))])
% disp(['RCQuEST: ', num2str(Prial_vec(10))])

end
