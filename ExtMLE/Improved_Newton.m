function [b_k tem d_k]= Improved_Newton(U,table)
% Newton method for Poisson
% Numerical problems: don't use gb_k and gb_kk (see example_num_problem)
% input:
%       U       the full rank design matix
%       table   the observe table of counts

alpha=1e-4; % slope modification parameter
epsilon=1e-22;

lambda=3*epsilon; 
iter=0;

% Initialize
m_k=table+1;
b_k = inv(U'*U)*U' *log(m_k);

m_kk=m_k;
b_kk=b_k;

mu_k = log(m_k);
mu_kk=mu_k;

fprintf('Iteration\t   f(x)   \t    Newton Step norm \t  Line-Search parameter\n');
fprintf('---------\t----------\t   ------------------\t  ---------------------\n');


% Newton Steps
while((table' * mu_kk - sum(exp(mu_kk)) - table' * mu_k + sum(exp(mu_k)) > 0)||(iter==0))    
    % updates
    b_k=b_kk;
    mu_k=U*b_k;
    m_k = exp(mu_k);
    grad_k = U' * (table - m_k);
    
    % Negative of the hessian    
    nhess_k = Compute_Hessian(U,m_k);
    
    %Terminating criterion for Newton's method (see Boyd pag. 474)
    lambda=grad_k' * inv(nhess_k) * grad_k;
    iter=iter+1;
    if(lambda/2 <= epsilon)
        b_k=b_kk;
        fprintf('\nConvergence achieved with lambda^2/2: %f\n',lambda/2);
        return;
    end

    if(iter>=100)
        b_k = b_kk;
        fprintf('\nExceed number of iterations: %d\n',iter);
        return;
    end
% Compute Newton direction. 
% Here can do better using Cholesky decomposition since nhess_k is 
% positive definite.

    d_k = nhess_k \ grad_k;
    
% Find preferred direction using backtracking line search
% Uses Numerical Recipes for the lines searches and backtracking
% methodology
% gb_k is the value of the funciton at b_k
% dgd = gradient of the function times d_k (directional derivatiove of the 
% function at b_k in the direction d_k)
% See page 390 of:
% Numerical Recipes in C++: The Art of Scientific Computing 
% by William H. Press (Editor), Saul A. Teukolsky (Editor), William T.
% Vetterling, Brian P. Flannery 
% Cambridge University Press
% 2002

    f2=0;
    alam2=0;
    % scale d_k if attemped step is too big (not implemented)
    % .
    % .
    % .
    
    % Slope:
    dgd = grad_k' * d_k; 
    % If slope dgd is zero or negative exit
    if(dgd <= 0)
        error('Round-off problem. Exiting.')
    end
    
    % Compute alamin
    alamin = 1e-22;
    %
    alam=1;
    
    
    while(alam>0)
        if(alam < alamin)
            fprintf('\nExiting because convergence achieved in line search: %f\n',alam);       
            tem=b_k;
            b_k = b_kk;
            b_kk = tem;
            return;
        end
        
        b_kk = b_k + alam .* d_k;
         mu_kk = U*b_kk;       
 
        % If sufficient function increase break out from inner loop
 
        if(table' * mu_kk - sum(exp(mu_kk)) - table' * mu_k + sum(exp(mu_k)) - (alpha * alam * dgd) >=0)
            table' * mu_kk - sum(exp(mu_kk)) - table' * mu_k + sum(exp(mu_k));
            break;
        end

        
        if(alam==1) % First time
            tmplam= - dgd/(2 * (table' * mu_kk - sum(exp(mu_kk)) - table' * mu_k + sum(exp(mu_k)) - dgd));          
        else

            rhs1 = ( table' * mu_kk - sum(exp(mu_kk)) - table' * mu_k + sum(exp(mu_k)) ) - alam * dgd;
            rhs2 = f2 - (  table' * mu_k - sum(exp(mu_k))  ) - alam2 * dgd;
            a = ( rhs1/(alam^2) - rhs2/(alam2^2) )/( alam - alam2 );
            b = (-alam2*rhs1/(alam^2) + alam*rhs2/(alam2^2) )/( alam - alam2 );          
            if(a==0) 
                tmplam = -dgd/(2*b);
            else
                disc=b^2 - 3*a*dgd;
                if(disc<0) 
                    tmplam = .5 * alam;
                elseif(b <=0)                 
                    tmplam = (-b + sqrt(disc))/(3*a);
                else                    
                    tmplam = -dgd/(b + sqrt(disc));
                end                
            end 
            
            if(tmplam > 0.5 * alam)                
                tmplam = 0.5 * alam;
            end
        end 
        
        alam2=alam;
        f2=table' * mu_kk - sum(exp(mu_kk));
        alam = max(tmplam, 0.1*alam);
        
    end 
    
    fprintf('  %d\t\t%f\t\t%f\t\t%f\n',iter,table' * mu_kk - sum(exp(mu_kk)),norm(alam .* d_k),alam);
    
    if(isequal(b_k,b_kk))
        disp('Convergence achieved in the parameter values...Exiting')
       return;
    end            

end


if( table' * mu_kk - sum(exp(mu_kk)) - table' * mu_k + sum(exp(mu_k)))
    disp('Function not increasing anymore...')
end




function hes = Compute_Hessian(A,m)
% Compute A' * diag(m) * A

[r,c]=size(A);
AA = zeros([r c]);
for i = 1:c
   AA(:,i) = A(:,i).*m;  
end
hes = AA' * A;