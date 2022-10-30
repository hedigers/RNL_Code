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
    Z=x./sqrt(sum(x.^2,2));
    [~,eig_fix] = QIS(Z,1);
else
    [~,eig_fix] = QIS(x,1);
    Z=x;
end

[eig_fix] = sortLambda(eig_fix, n, p);


if complicated==true
    V=VIteration(Z,eig_fix);
end

Hinv=V*diag(eig_fix.^(-1))*V';
lower = diag((Z*Hinv*Z')./p);
Y=(Z./sqrt(lower));
H0=QIS(Y,1);

Sigma=diag(sqrt(dsquare))*H0*diag(sqrt(dsquare));
Sigma = SymPDcovmatrix(Sigma);

Sigma=p.*Sigma./trace(Sigma);

end

