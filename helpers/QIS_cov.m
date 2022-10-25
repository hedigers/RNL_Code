
function [sigmahat,deltaQIS]=QIS_cov(sample,T,p) 

n=T-1;                              % adjust effective sample size
c=p/n; 
[u,lambda]=eig(sample,'vector');    % spectral decomposition
[lambda,isort]=sort(lambda);        % sort eigenvalues in ascending order
u=u(:,isort);                       % eigenvectors follow their eigenvalues
%%% COMPUTE Quadratic-Inverse Shrinkage estimator of the covariance matrix %%%
h=min(c^2,1/c^2)^0.35/p^0.35;             % smoothing parameter
invlambda=1./lambda(max(1,p-n+1):p);      % inverse of (non-null) eigenvalues
Lj=repmat(invlambda,[1 min(p,n)])';       % like  1/lambda_j
Lj_i=Lj-Lj';                          % like (1/lambda_j)-(1/lambda_i)
theta=mean(Lj.*Lj_i./(Lj_i.^2+h^2.*Lj.^2),2);     % smoothed Stein shrinker
Htheta=mean(Lj.*(h.*Lj)./(Lj_i.^2+h^2.*Lj.^2),2); % its conjugate
Atheta2=theta.^2+Htheta.^2;                      % its squared amplitude
if p<=n % case where sample covariance matrix is not singular
   delta=1./((1-c)^2*invlambda+2*c*(1-c)*invlambda.*theta ...
      +c^2*invlambda.*Atheta2);           % optimally shrunk eigenvalues
else % case where sample covariance matrix is singular
   delta0=1./((c-1)*mean(invlambda));     % shrinkage of null eigenvalues
   delta=[repmat(delta0,[p-n 1]);1./(invlambda.*Atheta2)];
end
deltaQIS=delta.*(sum(lambda)/sum(delta)); % preserve trace
% deltaQIS=delta;
sigmahat=u*diag(deltaQIS)*u'; 

if p > T
    eig_fix_0=deltaQIS(1:(p-T+1));
    un=unique(eig_fix_0);
    if length(un) > 1

        nr=zeros(length(un),1);
        for i=1:length(un)
            nr(i)=length(find(eig_fix_0==un(i)));
        end

        [~,ind]=max(nr);
        % Find actual value
        val=eig_fix_0(eig_fix_0==un(ind));
        eig_fix_0=ones(length(eig_fix_0),1)*val(1);

    end

    eig_fix_1=sort(deltaQIS((p-T+2):end)); %eig_fix((p-T+2):end); %
    deltaQIS=[eig_fix_0; eig_fix_1];

else

    deltaQIS = sort(deltaQIS);

end


end