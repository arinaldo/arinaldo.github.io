function coface = Compute_Facial_Newton(B,maxiter)
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
%   It just requires a modification of the routine Improved_Newton.m.
%   Rounding to zero for small elements of mu is performed in order to
%   compute the facial set.
% Inputs:
%   B: matrix sch that size(B,2) = number of parameters under Poisson
%                                  sampling.
% Outpus:
%   coface: index set containing the coface with maximal support.
%           coface is a subset of 1:size(B,1).
%   exitflag: 1 convergence achieved
%             0 exceeded number of iterations




if(nargin==1)
    maxiter = 100;    
end

muzero=1e-8;

alpha=1e-4; % slope modification parameter
epsilon=1e-22;
iter=1;

exitflag = 1;


% Line search parameters:

f2=0;
alam2=0;

% Initialize values for b_k
b_k = ones(size(B,2),1)./(size(B,1)*100);
mu_k = B*b_k;


mu_kk=mu_k;
m_k = exp(mu_k);
m_kk=m_k;
b_kk=b_k;


 fprintf('Iteration\t   f(x)   \t    Newton Step norm \t  Line-Search parameter\n');
 fprintf('---------\t----------\t   ------------------\t  ---------------------\n');



%======================  Beginning of Newton Steps  =======================
%lambda=3*epsilon; 
%while(lambda/2 > epsilon)
%while(gb_kk > gb_k)
%while( - sum(exp(mu_kk)) + sum(exp(mu_k)) > 0)
while( (- sum(exp(mu_kk)) + sum(exp(mu_k)) > 0) || (iter==1) )
       
    b_k=b_kk;
    mu_k=B*b_k;
    m_k = exp(mu_k);
    
    
    % ====  Return a coface, if any is found.  ==== 
    %if(norm(b_k - b_kk,inf)) < 1e-5
    if((norm(m_k - m_kk,inf) < 1e-8)&&(iter>1))
        % Make small values zeros to avoid rounding errors.    
        mu_temp=mu_k;
        mu_temp(abs(mu_k)<muzero)=0;
        if(all(mu_temp<=0)||all(mu_temp>=0))
            coface=find(abs(mu_temp)>0);
            return;
        else
            coface=[];
            return;
        end        
        
    end
    % =============================================
     
    
    %disp(['New step. gb: ', num2str(  - sum(exp(mu_k)))]);
    
    % Compute gradient. 
    grad_k = B' * (- m_k);

    % Compute Hessian.
    hess_k = B' * diag(m_k) * B;

    %Terminating criterion for Newton's method (see Boyd pag. 487)
    lambda=grad_k' * inv(hess_k) * grad_k;
    %disp([num2str(iter),': ',num2str(lambda/2)]);
    iter=iter+1;
    if(lambda/2 <= epsilon)        
        fprintf('\nConvergence achieved with lambda^2/2: %f\n',lambda/2);
        
        % Return a coface, if any is found.
        % Make small values zeros to avoid rounding errors.    
        mu_temp=mu_k;
        mu_temp(abs(mu_k)<muzero)=0;
        if(all(mu_temp<=0)||all(mu_temp>=0))
            coface=find(abs(mu_temp)>0);
            return;
        else
            coface=[];
            return;
        end  
    end

    % Too many iterations: exiting
    if(iter>=maxiter)
        fprintf('\nExceed number of iterations: %d\n',iter);
        exitflag = 0;
        
        % Return a coface, if any is found.
        % Make small values zeros to avoid rounding errors.    
        mu_temp=mu_k;
        mu_temp(abs(mu_k)<muzero)=0;
        if(all(mu_temp<=0)||all(mu_temp>=0))
            coface=find(abs(mu_temp)>0);
            return;
        else
            coface=[];
            return;
        end  
        
    end
    
    
    
    % Compute Newton direction. 
    % Here can do better using Cholesky decomposition since hess_k is 
    % positive definite.
    d_k = hess_k \ grad_k;

    
    
    % Find preferred direction using backtracking line search
    % methodology
    % gb_k is the value of the funciton at b_k
    % dgd = gradient of the function times d_k (directional derivative 
    % of the function at b_k in the direction d_k)
    % See Boyd page 464 (t there is alam here. alpha there is alpha here. 
    % beta there is not used here because I use instead the method 
    % on page 390 of Numerical Recipes in C++ for backtracing.

    
    % scale d_k if attemped step is too big (not implemented)
    % .
    % .
    % .
    
    
    % Slope:
    dgd = grad_k' * d_k; 
    % If slope dgd is zero or negative exit
    if(dgd <= 0)
        error('Round-off problem in Imporved_Newton. Exiting.')
        coface=-1;
        exitflag = -1;
    end
    
    % Compute alamin
    alamin = 1e-22;
    
    alam=1;
    
    %====== Newton step with backtracking using cubic approximation. ======
    while(alam>0)
        %disp('alam :')
        %disp(alam)
        % Convergence along d_k achieved (to be verified in outer loop)
        if(alam < alamin)                        
            fprintf('\nExiting because convergence achieved in line search: %f\n',alam);
         %   disp(isequal(b_k, b_kk))
            tem=b_k;
            b_k = b_kk;
            b_kk = tem;
            
            
            % Return a coface, if any is found.
            % Make small values zeros to avoid rounding errors.    
            mu_temp=mu_k;
            mu_temp(abs(mu_k)<muzero)=0;
            if(all(mu_temp<=0)||all(mu_temp>=0))
                coface=find(abs(mu_temp)>0);
                return;
            else
                coface=[];
                return;
            end  
            
            
        end
        
       
        % New b_kk
        b_kk = b_k + alam .* d_k;
        % value of g at b_kk
        %gb_kk = table' * (mu_kk) - sum(exp(mu_kk));
        mu_kk = B*b_kk;       
        
    
        % If sufficient function increase break out from  while(alam>0) 
        % loop. Try full Newton step first.
        % if(gb_kk >= gb_k + (alpha * alam * dgd))
        if(- sum(exp(mu_kk)) + sum(exp(mu_k)) - (alpha * alam * dgd) >=0)
            %disp(['Increase in gb: '])%,num2str(gb_k),' and gb_kk: ', num2str(gb_kk)])                        
            %disp([' alam: ', num2str(alam)]);
            %disp(['Norm increase: ',num2str(norm(alam .* d_k)), ' Incremenetal variation of gb_kk: ',num2str(alpha * alam * dgd)]);
            %disp((alpha * alam * dgd)>0);
            %disp(' ' );
            %disp('');
            break;
        end
        
        % Backtracking
        %disp('Changing alam')
        
        if(alam==1) % First time
            %tmplam= - dgd/(2 * (gb_kk - gb_k - dgd));          
            tmplam= - dgd/(2 * (- sum(exp(mu_kk)) + sum(exp(mu_k)) - dgd));          
        else
            %rhs1 = gb_kk - gb_k - alam * dgd;
            rhs1 = ( - sum(exp(mu_kk)) + sum(exp(mu_k)) ) - alam * dgd;
            
            %rhs2 = f2 - gb_k - alam2 * dgd;
            rhs2 = f2 - (  - sum(exp(mu_k))  ) - alam2 * dgd;
            
            
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
            end % close if(a==0)            
            
            if(tmplam > 0.5 * alam)                
                tmplam = 0.5 * alam;
            end
        end % close if(alam==1)
        
        
        
        alam2=alam;
        %f2=gb_kk;
        f2= - sum(exp(mu_kk));
        alam = max(tmplam, 0.1*alam);
        
    end % close while(alam>0)
    %======================================================================
    fprintf('  %d\t\t%f\t\t%f\t\t%f\n',iter,- sum(exp(mu_kk)),norm(alam .* d_k),alam);
    if(isequal(b_k,b_kk))
        disp('Convergence achieved in the parameter values...Exiting')
        
        % Return a coface, if any is found.
        % Make small values zeros to avoid rounding errors.    
        mu_temp=mu_k;
        mu_temp(abs(mu_k)<muzero)=0;
        if(all(mu_temp<=0)||all(mu_temp>=0))
            coface=find(abs(mu_temp)>0);
            return;
        else
            coface=[];
            return;
        end  
        
    end
  
end


if( - sum(exp(mu_kk)) + sum(exp(mu_k)))
    disp('Function not increasing anymore...')
end


% Return a coface, if any is found. TO DO  replace 0 with myzero
if(all(mu_k<=0)||all(mu_k>=0))
    coface=find(abs(mu_k)>0);
else
    coface=[];
end

