
% Percentage Relative Improvement in Average Loss (PRIAL)

%% Input:

% S: sample covariance matrix

% Sigma_hat: estimated covariance matrix

% Sigma: true covariance matrix

% note: the provided matrices can also be 3 dimensional arrays of multiple simulation
% runs.

%% Output:

% prial: numerical value of the prial.

%%

function [prial] = PRIAL(S,Sigma_hat,Sigma)

for i = 1:size(Sigma_hat,3)

    if size(Sigma,3) == 1
        a(i) = norm(Sigma_hat(:,:,i)-Sigma,'fro')^2;
        b(i) = norm(S(:,:,i)-Sigma,'fro')^2;
    else
        a(i) = norm(Sigma_hat(:,:,i)-Sigma(:,:,i),'fro')^2;
        b(i) = norm(S(:,:,i)-Sigma(:,:,i),'fro')^2;
    end
    
end

prial = (1-mean(a)/mean(b))*100;


