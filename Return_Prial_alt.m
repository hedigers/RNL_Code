function [Prial_vec] = Return_Prial_alt(p, n, mu, H, lam, S)

rng(747, 'twister');

alpha_vec = 0:0.2:1;

H = SymPDcovmatrix(H);
H_true = H/trace(H)*p;

H_sample_s = NaN(p,p,S);
H_rob_nl_s = NaN(p,p,S,length(alpha_vec)+2);

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

    % iterating through alphas
    for i = 1:length(alpha_vec)
        H_rob_nl = RNL(X, true, true, false, alpha_vec(i),false);
        H_rob_nl = H_rob_nl/trace(H_rob_nl)*p;
        H_rob_nl_s(:,:,s,i) = H_rob_nl;
    end

    % R-NL
    H_rob_nl = RNL(X, false, true, true, 1, true);
    H_rob_nl = H_rob_nl/trace(H_rob_nl)*p;
    H_rob_nl_s(:,:,s,i+1) = H_rob_nl;

    % only V1
    H_rob_nl = RNL(X, false, false, true, 1, true);
    H_rob_nl = H_rob_nl/trace(H_rob_nl)*p;
    H_rob_nl_s(:,:,s,i+2) = H_rob_nl;

end

% Standard PRIAL
Prial_vec = NaN(length(alpha_vec)+2,1);
for i = 1:(length(alpha_vec)+2)
    Prial_vec(i) = PRIAL(H_sample_s,H_rob_nl_s(:,:,:,i),H_true);
end

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
