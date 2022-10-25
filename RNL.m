function [Sigma] = RNL(X, complicated, normscale, corrbased)

if nargin < 3
   complicated=false; 
end

[n,p] = size(X);

location=mean(X);
x=X-location;

dsquare=var(x);
if corrbased == 0
    dsquare = ones(1,p);
end

x=x./sqrt(dsquare);
if normscale == 1
    z=x./sqrt(sum(x.^2,2));
    [~,eig_fix] = QIS(z,1);
else
    [~,eig_fix] = QIS(x,1);
    z=x;
end

eig_fix=sortLambda(eig_fix, n, p);

Psi_tmp = cov(z);
[V0,lam] = eig(Psi_tmp,'vector');
[~,ind] = sort(lam);
V0 = V0(:,ind);

if complicated==true
    Gamma=VIteration(z,eig_fix, V0);
end

Sigma=diag(sqrt(dsquare))*Gamma*diag(sqrt(dsquare));
Sigma = SymPDcovmatrix(Sigma);

end

