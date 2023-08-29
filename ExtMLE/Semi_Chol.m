function [R P ran]= Semi_Chol(A,piv,tole)
% Semi_Chol     Cholesky decomposition of a symmetric semi-definite via pivoting
%       Inputs  A   : symmetric semi-definite matrix
%               piv : 0-1 flag. if piv=0 (default), no pivoting is
%                     performed, in which case the function might return an error
%                     if the A is not definite positive
%               tole: tolerance for small values of the sqaured digonal elements of
%                     the Cholesky score
%       Output  R :   Upper triangular matrix
%               P :   Permutation vector. If length(P)=ran<size(A,1) then, A is
%                     semi-definite and ran is the rank of A.
%                     In addition, if A=U_0' * U_0, where U_0 is the U basis with removed 
%                     facial zeroes, then P(1:ran) contains the columns of U_0 spanning the
%                     restricted model.
%       Note:  A(P,P) - R'*R = 0 (the null matrix) or, more exactly, because of rounding errors, 
%              just a small one.

n= size(A,1);
P= 1:n;
R=zeros(n, n);
AA=A;
diagonal  = sub2ind([n,n], 1:n, 1:n); 

if(nargin<2)
    piv=0;
end

if(nargin==3)
    tol=tole;
else
    tol = 1e-12;
end

for k=1:n
    if(piv==1)
        
        % Get biggest element in AA diagonal
        [a_max, k_max]=max(AA(diagonal(k:n)));
        
        if(k_max+k-1 ~= k)
            % pivoting
            % update the permutations vector
            temp=P(k);
            P(k)=P(k_max+k-1);
            P(k_max+k-1)=temp;
            % permute appropriate rows and columns of AA and R (can be
            % improved)
            AA([k (k_max+k-1)],:)=AA([(k_max+k-1) k],:);
            AA(:,[k (k_max+k-1)])=AA(:,[(k_max+k-1) k]);
            R([k (k_max+k-1)],:)=R([(k_max+k-1) k],:);
            R(:,[k (k_max+k-1)])=R(:,[(k_max+k-1) k]);        
            
			% ALTERNATIVE PERMUTING
            %AA([k (k_max+k-1)],k:n)=AA([(k_max+k-1) k],k:n);
            %AA(k:n,[k (k_max+k-1)])=AA(k:n,[(k_max+k-1) k]);
            %R([k (k_max+k-1)],k:n)=R([(k_max+k-1) k],k:n);
            %R(k:n,[k (k_max+k-1)])=R(k:n,[(k_max+k-1) k]);
            
        end
           
    end
    % Proceed with usual Cholesky decomposition
    v=sqrt(AA(k,k));
    if(piv==1)
        % If max element in diagonal AA is less than tol end
        if(v <= tol)
            P=P(1:k-1);
            break;
        end
    end
    R(k,k)=v;
    R(k,k+1:n)=AA(k,k+1:n)/R(k,k);
    AA(k+1:n,k+1:n)=AA(k+1:n,k+1:n) - R(k,k+1:n)'*R(k,k+1:n);
    
end
ran=length(P);
P=[P sort(setdiff((1:n),P))];