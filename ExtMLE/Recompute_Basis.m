function X = Recompute_Basis(U,chol)
% Recompute_Basis: given a matrix U not of full-column rank, it produces a
%   matrix of full columns rank whose columns are a subset of the columns
%   of U.
% Inputs:
%   U: non- full rank design matrix.
%   chol: if different than 0 enforces the usage of the Cholesky
%         decomposition. Otherwise, rref is used. Deafualt = 0.
% OutputS:
%   X: full column-rank matrix whose columns form a subset of the columns 
%      of the original matrix U.

if (nargin<2)
    chol=0;
end

A=U'*U;
if(~chol)
    % Use rref.
    [R jb]=rref(A);
    X=U(:,jb);
else
    % Use Cholesky decomposition.
    [R P ran]=Semi_Chol(A); 
    X=U(:,P(1:ran));
end