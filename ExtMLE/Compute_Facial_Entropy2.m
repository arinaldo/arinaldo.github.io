function [coface,exitflag]= Compute_Facial_Entropy2(B,tol)
%function [coface, exitflag, output] = Compute_Facial_Entropy2(B,tol)
% Compute_Facial_Entropy1: compute the facial set using Entropy.
%   B is a integer matrix such that the MLE exists if and only if the
%   system Bx \gneq 0 has a solution. The co-face corresponds to a
%   solution with maximal support.
%   This nonlinear program solves the problem:
%
%   min H(x)
%   s.t. B'x =  0
%        x   >= 0
%        1'x =  1
%   
%   where B = Uz * X, with X being a matrix whose columns (integer vectors) 
%   form a basis for kernel(Uznot). Zero rows from B are already expunged.
%   H(x) is the entropy function.
%   If x* is the optimal solution then, zero entries in supp(Bx*) give the 
%   coordinates for the coface.
%   Rounding to zero of small elements of Bx is performed in 
%   order to compute the facial set.
% Inputs:
%   B: matrix such that size(B,2) = number of parameters under Poisson
%                                   sampling.
%   tol: tolerance level for the rounding (default is 1e-8).
% Outputs:
%   coface: index set containing the coface with maximal support.
%           coface is a subset of 1:size(B,1).
%
% Comments: the following options are used for linprog:
%           optimset('Simplex','off','LargeScale','on','Display','iter')
% so that interior points method is used instead of the simplex method.



if(nargin<2)
    tol=1e-12;
end


% Set the optimization options for optimset.
options = optimset('GradObj','on','Hessian','on','Display','iter','FunValCheck','on','Diagnostics','on');

% Setting up parameters for the equalitites and inequalities consraints.
rB=size(B,1);
cB=size(B,2);       

Aeq = [B';ones(1,rB)];
beq=[zeros(cB,1);1];
lb = zeros(rB,1);
% This constraint could be incorporated as lb = zeros(rB,1);
%A=[-eye(rB)];
%b=zeros(rB,1);


% Initial value for x (infeasible).
x0=1/rB .* ones(rB,1);

% Use fmincon to find the optimal value out.
%[out, Entropy2_val, exitflag, output] = fmincon(@Entropy2, x0, A, b, Aeq, beq,[],[],[],options);
% Using constraint lb = zeros(rB,1), the call would be:
[out, Entropy2_val, exitflag, output] = fmincon(@Entropy2, x0, [], [], Aeq, beq,lb,[],[],options);


if((exitflag==-2)||(coface==0))
    coface = 1:rB;
else
    coface = find(out <= tol);        
end

% Nested function for the Entropy.
    function [f,Grad,Hess] = Entropy2(x)            
        f=0;    
        if(all(x<=0))
            f=0;
            if(nargout>1)     
                Grad=Inf .* ones(size(x));
                if(nargout>2)
                    Hess = zeros(length(x));
                end
            end            
        else
            index  = find(x>0);       
            f = f+ x(index)' *log(x(index));
            if(nargout>1)            
                % Gradient as a column vector.
                Grad = zeros(size(x));
                Grad(index) = 1+log(x(index));
                if(nargout>2)
                    % Hessian.
                    Hess=zeros(size(x));
                    Hess(index) = 1./x(index);
                    Hess = diag(Hess);
                end
            end
        end
    end


end