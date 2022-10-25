function A=SymPDcovmatrix(A,tol)
% function corrects a covariance matrix A to be symmetric positive definite 
% it uses eigenvalue decomposition and shifts all small eigenvalues to tol.
% INPUTS: 
% A - matrix
% tol - (optional, default tol = 1e-04) minimum value for all eigenvalues
% in the output
% OUTPUTS: A - corrected matrix A.

[n,m]=size(A);% make sure it is a square matrix
if ~(n==m)
    error('Input matrix has to be a square matrix '), end
if nargin<2
    tol=1e-04;
end
% symmetrize the matrix.
%A=A-diag(diag(A))+eye(n,m);
A=(A+A')/2;

% correct it eigenvalues to make sure the matrix is positive definite.
[V,D]=eig(A); seig=diag(D); bad=find(seig<tol);
if ~isempty(bad), seig(bad)=tol; D=diag(seig); A=V*D*V';end


