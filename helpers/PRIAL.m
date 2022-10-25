function [prial] = PRIAL(S,Sigma_hat,Sigma)

p = size(Sigma,1);

for i = 1:size(Sigma_hat,3)

    if size(Sigma,3) == 1
        a(i) = norm(Sigma_hat(:,:,i)-Sigma,'fro')^2;
        b(i) = norm(S(:,:,i)-Sigma,'fro')^2;
    else
        a(i) = norm(Sigma_hat(:,:,i)-Sigma(:,:,i),'fro')^2;
        b(i) = norm(S(:,:,i)-Sigma(:,:,i),'fro')^2;
    end

    %  a(i) = norm(Sigma_hat(:,:,i)-Sigma)^2;
    %  b(i) = norm(S(:,:,i)-Sigma)^2;
    %     a(i) = (norm(Sigma_hat(:,:,i)-Sigma,'fro')/sqrt(p))^2;
    %     b(i) = (norm(S(:,:,i)-Sigma,'fro')/sqrt(p))^2;
    %
    %     A = Sigma_hat(:,:,i)-Sigma;
    %     B = S(:,:,i)-Sigma;
    %     a(i) = sqrt(trace(A'*A)/sqrt(p))^2;
    %     b(i) = sqrt(trace(B'*B)/sqrt(p))^2;
end

prial = (1 - mean(a)/mean(b)) * 100;


