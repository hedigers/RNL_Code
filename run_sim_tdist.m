clear

rng(747, 'twister');

n = 100; % choose the number of observations

S = 1000;
lam_vec = [3:15,30,60,90,120,240,500,inf];
p_vec = [200,400]; % choose the number of dimensions
II = length(p_vec);
HH = 6;
Prial_cell = cell(II,1);

Prial_mat = NaN(length(lam_vec),HH,II,8);
% Note: the last index stands for the number of different estimators.

% uncomment the below to use a parpool structure
% parpool(length(lam_vec)+1, 'IdleTimeout', Inf)

for ii = 1:II

    p = p_vec(ii);
    mu = zeros(p,1);

    for hh = 1:HH

        if hh == 1

            %%%%%% Identity %%%%%%%%%
            H = eye(p);

        elseif hh == 2

            %%%%%%% AR %%%%%%%%%%%%%%
            r = 0.7;
            for i = 1:p
                for j = 1:p
                    H(i,j) = r^abs(i-j);
                end
            end

        elseif hh == 3

            %%%%%%% Full Matrix %%%%%%%%%%%%%%%
            H = 0.5*ones(p,p);
            for i = 1:p
                H(i,i) = 1;
            end

        elseif hh == 4

            %%%%%%% Base %%%%%%%%%%%%%%%
            if p == 3
                H = diag([10,3,1]);
            else
                H = diag([repmat(10,1,0.4*p) repmat(3,1,0.4*p) repmat(1,1,0.2*p)]);
            end

        elseif hh == 5

            %%%%%%% AR (non-constant diag) %%%%%%%%%%%%%%%
            r = 0.7;
            for i = 1:p
                for j = 1:p
                    H(i,j) = r^abs(i-j);
                end
            end
            tmp = [repmat(10,1,0.4*p) repmat(3,1,0.4*p) repmat(1,1,0.2*p)];
            H = diag(sqrt(tmp))*H*diag(sqrt(tmp));
            H = nearestSPD(H);

        elseif hh == 6

            %%%%%%%%%%%% Full Matrix (non-constant diag) %%%%%%%%%%%%%%%%%%
            H = 0.5*ones(p,p);
            for i = 1:p
                H(i,i) = 1;
            end

            tmp = [repmat(10,1,0.4*p) repmat(3,1,0.4*p) repmat(1,1,0.2*p)];
            H = diag(sqrt(tmp))*H*diag(sqrt(tmp));
            H = nearestSPD(H);

        end

        disp(['hh = ', num2str(hh)])

       for idx = 1:length(lam_vec) % you can use parfor here

            lam = lam_vec(idx);
            Prial_vec = Return_Prial(p, n, mu, H, lam, S);
            Prial_mat(idx,hh,ii,:) = Prial_vec;

        end

        save(['Prial_mat_',num2str(n),'_',num2str(S),'.mat'], 'Prial_mat','-v7.3')

    end

end
