function [beta, exitflag, output] = Matlab_Newton(U,table)
%[beta, exitflag, output] = Matlab_Newton(U,table)
% Maximize the Poisson Log-Likelihood function using the
% Matlab function fminunc.
% Input: U      design matrix
%        table  observed counts
% Output: beta      MLE of the parameter
%         exitflag  flags for fminunc function
%         output    info from the optimization from fminunc


% Set the Gradj and Hessian options on.
options = optimset('GradObj','on','Hessian','on','Display','iter');

% Initial value for beta
beta0=ones(size(U,2),1);
% beta0 = inv(U'*U)*U'*table;

[beta, LogLik_val, exitflag, output] = fminunc(@LogLik_Pois, beta0, options);


% Nested function for the Log-Likelihood.
% It uses variables U and table.
    function [f,Grad,Hess] = LogLik_Pois(b)            
        mu=U*b;
        f = - table' * mu + sum(exp(mu));
        if(nargout>1)            
            % Gradient
            Grad = U' * (exp(mu) - table);
            if(nargout>2)
                % Hessian
                %Hess = Compute_Hessian(U,exp(mu)); %U' * diag(exp(mu)) * U;
                [r,c]=size(U);
                UU = zeros([r c]);
                for i = 1:c
                    UU(:,i) = U(:,i).*exp(mu);  
                end
                Hess = UU' * U;
            end
        end        
    end




end