function flag = Has_MLE(counts,Table,Faces)
% Has_MLE: decide whether it has a MLE or ,not.
%   The function implements the following LP procedure in 
%   1 + length(counts) variables.
%
%   min -s
%   s.t. s >=0
%        A x = b
%        x_i - s >= 0
%
%   The MLE exists if and only if the optimal solution is such that s > 0.
%   Rounding to zero are applied.
% Inputs:
%   counts: table counts.
%   Table: structure encoding the Log-Linear Model.
%   Faces: class of subsets holding the interactions (optional).
% Output:
%   flag: 1 is the MLE exists and 0 otherwise.


flag = 1;

if(nargin<2)
   error('Function Has_MLE requires two inputs. Exiting...'); 
end


if(nargin<3)   
    U = Make_Basis(Table);
else      
    U = Make_Basis(Table,'C',Faces);
end

% Setting up the LP parameters.

options=optimset('Simplex','off','LargeScale','on','Display','iter');
tol = 1e-8;

rU = size(U,1);

Aeq = [U' ones(size(U,2),1)];
beq = U' * counts;
A = [-eye(rU) ones(rU,1); [zeros(1,rU) -1] ];
b = zeros([rU+1 1]);
f = [zeros(rU,1);-1];

[out,fval,exitflag,output]=linprog(f,A,b,Aeq,beq,[],[],[],options);       

out(end)

if(out(end)<= tol)
    flag = 0;
end
