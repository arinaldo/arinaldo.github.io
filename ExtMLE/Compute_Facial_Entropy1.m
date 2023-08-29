function [coface,exitflag, output,out] = Compute_Facial_Entropy1(B,tol)
% Compute_Facial_Entropy1: compute the facial set using Entropy.
%   B is a integer matrix such that the MLE exists if and only if the
%   system Bx \gneq 0 has a solution. The co-face corresponds to a
%   solution with maximal support.
%   This nonlinear program solves the problem:
%
%   min H(x)
%   s.t. 1'Bx =  1
%        Bx   >= 0
%   
%   where B = Uz * X, with X being a matrix whose columns (integer vectors) 
%   form a basis for kernel(Uznot). Zero rows from B are already expunged.
%   H(x) is the entropy function.
%   If x* is the optimal solution then, supp(Bx*) gives the coordinates for
%   the coface.
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




% Setting up default value for the tolerance.
if(nargin<2)
    tol=1e-8;
end


% Set the Gradj and Hessian options on.ptimset
options = optimset('GradObj','on','Hessian','on','Display','iter','FunValCheck','on','Diagnostics','on');

% Setting up parameters for the equalitites and inequalities consraints.
rB=size(B,1);
cB=size(B,2);       

%Aeq = ones(1,rB)*B;
%beq=1;
%A=-B;
%b=zeros(rB,1);

d = ones(1,rB)*B;
A=[-B;d;-d];
b = [zeros(rB,1);1;-1];


% Initial value for x (infeasible).
x0=ones(cB,1);

%[out, Entropy1_val, exitflag, output] = fmincon(@Entropy1, x0, A, b, Aeq, beq,[],[],[],options);
[out, Entropy1_val, exitflag, output] = fmincon(@Entropy1, x0, A, b,[],[],[],[],[],options);

%y = B*out;
if((exitflag ==-2)||(exitflag==0))
    coface = [];
    return;
else
    coface = find(B*out>tol);
end

% Nested function for the Entropy..
    function [f,Grad,Hess] = Entropy1(x)            
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
            xx=B*x;
            index=find(xx>0);
            f = f + xx(index)' *log(xx(index));
            if(nargout>1)            
                % Gradient (column vector).            
                Gradxx = zeros(size(xx));
                Gradxx(index) = 1+log(xx(index));
                Grad = B' * Gradxx;
                if(nargout>2)
                    % Hessian.
                    Hessxx=zeros(size(xx));
                    Hessxx(index) = 1./xx(index);
                    Hess = B' * diag(Hessxx) * B;
                end
            end
        end
    end


end
