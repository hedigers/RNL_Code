function [eig_fix] = sortLambda(eig_fix, T, N)

%%%% Note: The first N-T elements in eig_fix necessarily correspond to the
%%%% zero eigenvalues.
eig_fix = sort(eig_fix);
if N > T
    % First eigenvalues=Smallest eigenvalues (in this case corresponding
    % to the sample eigenvalues with zeros)
    eig_fix_0=eig_fix(1:(N-T+1));
    un=unique(eig_fix_0);
    if length(un) > 1

        nr=zeros(length(un),1);
        for i=1:length(un)
            nr(i)=length(find(eig_fix_0==un(i)));
        end

        [~,ind]=max(nr);
        eig_fix_0=ones(length(eig_fix_0),1)*un(ind);

    end

    eig_fix_1=sort(eig_fix((N-T+2):end));
    eig_fix=[eig_fix_0; eig_fix_1];


end
end