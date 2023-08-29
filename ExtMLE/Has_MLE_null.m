function flag = Has_MLE_null(counts,B0)
% Has_MLE: decide whether it has a MLE or ,not.
%   The function implements the following LP procedure in 
%   1 + size(B0,2) variables.
%
%   min -s
%   s.t. s >=0
%        B0 x - s 1 >= 0
%       
%
%   The MLE exists if and only if the optimal solution is such that s > 0.
%   Rounding to zero are applied.
% Inputs:
%   counts: table counts.
%   B0: sub-matrix of B = kernel(U') whose rows corresponds to the zero 
%       cells.
% Output:
%   flag: 1 is the MLE exists and 0 otherwise.


flag = 1;

if(nargin<2)
    error('Function Has_MLE requires two inputs. Exiting...'); 
end


% Setting up the LP parameters.

options=optimset('Simplex','off','LargeScale','on','Display','iter');

tol = 1e-8;
[rB0, cB0] = size(B0);

b = zeros(rB0+1,1);
A = [-B0 ones(rB0,1); zeros(1,cB0) -1 ];
f = [zeros(cB0,1);-1];

[out,fval,exitflag,output]=linprog(f,A,b,[],[],[],[],[],options);       

out(end)

if(out(end)<= tol)
    flag = 0;
end