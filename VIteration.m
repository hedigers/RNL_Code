function [V] = VIteration(Z, eig_fix)

tol = 1e-4;
[~,p] = size(Z);
eig_fix_inv=eig_fix.^(-1);
%V=V0;
V=eye(p);

i=0;
diagnorm=1000;

while abs(diagnorm)> tol

    i=i+1;

%     if i > 200
%         warning("More than 200 iterations")
%         break
%     end

    V_old=V;

    Hinvold=V_old*diag(eig_fix_inv)*V_old';
    H_old=V_old*diag(eig_fix)*V_old';

    lower = diag((Z*Hinvold*Z')./p);
    ratioV = (Z./lower)'*(Z);
    ratioV=(ratioV+ratioV')/2;

    [V,lam]=eig(ratioV,'vector');
    [~, ind] =sort(lam);
    V=V(:,ind);

    diagnorm=norm(V'*ratioV*V*diag(eig_fix_inv) -diag(eig_fix_inv)*V'*ratioV*V );
    %%%norm(H_old-V*diag(eig_fix)*V');

end



end