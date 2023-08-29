function beta = Compute_Facial_Matlab_Newton(B)
% Compute_Facial_Newton: compute the facial set using Newton-Rapson's
%   method.
%   B is a integer matrix such that the MLE exists if and only if the
%   system y=Bx \gneq 0 has a solution. The co-face corresponds to a
%   solution with maximal support.
%   This program solves the problem:
%
%   max -sum_i \exp{(b_i, x)}
%
%   where b_i is the i-th row of the matrix B and B = Uz * X, with X being
%   a matrix whose columns form a basis for kernel(Uznot). 
%   This is an unconstrained optimizaiton problem which can be solved using
%   Newton-Rapson's method. The gradient and Hessian are easy to compute.
%   This version ises the Matlab unconstrained minimzation function
%   fminunc.
% Inputs:
%   B: matrix sch that size(B,2) = number of parameters under Poisson
%                                  sampling.
% Outpus:
%   coface: index set containing the coface with maximal support.
%           coface is a subset of 1:size(B,1).

% TO DO: SET UP THE TOLRANCE LEVEL FOR THE CONVERGENCE


% Set the Gradj and Hessian options on.ptimset
options = optimset('GradObj','on','Hessian','on','Display','iter','OutputFcn',@outfun);

tolm = 1e-8;
muzero = 1e-8;

% Initial value for x.
x0=ones(size(B,2),1);

[beta, LogLik_val, exitflag, output] = fminunc(@Facial, x0, options);

mu = B*beta;
mu(abs(mu)<muzero)=0;
if(all(mu<=0)||all(mu>=0))    
    coface=find(abs(mu_temp)>0);
    return;
else
    coface=[];
    return;
end     

% Nested functiond.
% It uses variable B.
    function [f,Grad,Hess] = Facial(x)            
        y=B*x;
        f = sum(exp(y));
        if(nargout>1)            
            % Gradient
            Grad = B' * exp(y);
            if(nargout>2)
                % Hessian
                Hess = B' * diag(exp(y)) * B;
            end
        end
    end


% Output function
    function stop = outfun(x,optimvalues, state)        
        persistent xi_1;
        stop = false;
        switch state
            case 'init'
                xi_1 = x;                
            case 'iter'
                mui_1 = B * xi_1;
                mui = B *x;
                if(norm(exp(mui_1) - exp(mui),inf) < tolm)
                    stop = true;
                else
                   xi_1 = x;
                end
            otherwise                                
        end         
    end