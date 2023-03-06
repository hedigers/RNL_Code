
% Helper function used in Robust_Lin_Shrink, where X is nxp and assumed to
% be centered.

%%

function [C] = C_fixpoint(rho, X, C_0, sp)

if ~exist('sp')
    sp = 0;
end

tol = 0.0001;
[n,p] = size(X);

C = C_0;
C_old = 100*C_0;

while norm(C-C_old) > tol

    C_old = C;
    C_inv = C\eye(p);

    lower = diag((X*C_inv*X')./p);
    ratio = (X./lower)'*(X);

    if sp == 1
        C = (1-rho)*(ratio/n) + rho*eye(p);
    elseif sp == 0
        C = (1-rho)*(ratio/n) + rho*eye(p);
        C = p*(C/trace(C));
    elseif sp == 2

        ratio = zeros(p,p);
        for t = 1:n
            ratio = ratio + p/n*((1-rho)*X(t,:)'*X(t,:) + rho*norm(X(t,:))^2.*eye(p)./p)./...
                (((1-rho)*X(t,:)*C_inv*X(t,:)' + rho*norm(X(t,:))^2*trace(C_inv)/p));
        end

        C = p*(ratio/trace(ratio));
    end

end

end