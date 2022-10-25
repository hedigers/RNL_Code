function [sig] = sig_sc(rho, X)

% X is nXp
% X is assumed to be centered

eps = 0.01;

[n,p] = size(X);
c_N = p/n;

C_rho = C_fixpoint(rho,X,eye(p));
C_rho_inv = inv(C_rho);

squared_norm = diag(X*X')';

B = squared_norm > eps;

gamma_hat = sum(diag(X(B,:)*C_rho_inv*X(B,:)')./squared_norm(B)');

gamma_hat = gamma_hat/(sum(B)*(1-(1-rho)*c_N));

sig = (gamma_hat/((1-rho)-((1-rho)^2)*c_N)) * ((ones(1,p)*C_rho_inv*(C_rho-(rho*eye(p)))*C_rho_inv*ones(p,1))/((ones(1,p)*C_rho_inv*ones(p,1))^2));

end