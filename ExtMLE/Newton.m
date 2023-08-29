function [b_k tem grad_k]= Newton(U,table)
% Newton method for Poisson

alpha=.4; % slope modification parameter
beta=.5; % backtracking line search parameter
epsilon=1e-22;
eps=1e-22;
lambda=3*epsilon; 
iter=1;

% Initialize
m_k=table+1;

b_k = inv(U'*U)*U' *log(m_k);


%%%% BEGIN 

% Newton-Steps

while(lambda/2 > epsilon)
    
    
    grad_k = U' * (table - m_k);

    hess_k = U' * diag(m_k) * U;

    %Terminating criterion for Newton's method (see Boyd pag. 474)
    lambda=grad_k' * (hess_k \ grad_k);%inv(hess_k) * grad_k;
    %disp([num2str(iter),': ',num2str(lambda/2)]);
    iter=iter+1;
    if(lambda/2 <= epsilon)
        disp(['Exiting with lambda^2/2: ',num2str(lambda/2)]);
        break;
    end

    if(iter>=100)
        disp(['Exceed number of iterations: ',num2str(iter)]);
        break;
    end
% Compute Newton direction. 
% Here can do better using Cholesky decomposition since hess_k is 
% positive definite.

    d_k = hess_k \ grad_k;

% Find preferred direction using backtracking line search
% gb_k is the value of the funciton at b_k
% dgd = gradient of the function times d_k (directional derivatiove of the 
% function at b_k in the direction d_k)
% see Boyd page 450

    mu_k = U*b_k;
    %gb_k = table' * (mu_k) - sum(exp(mu_k));
    dgd = grad_k' * d_k; 
% Initialize parameters for the backtrack
    b_kk = b_k+d_k;
    alpha_dgd = alpha .* dgd;
    mu_kk = U*b_kk;
    %gb_kk = table' * (mu_kk) - sum(exp(mu_kk));
% Baktrack algorithm
  %  while(gb_kk < gb_k + alpha_dgd)        
    while(table' * (mu_kk) - sum(exp(mu_kk)) < table' * (mu_k) - sum(exp(mu_k)) + alpha_dgd)        
        d_k = d_k .* beta;
        alpha_dgd = alpha_dgd .* beta;
        b_kk = b_k + d_k;
        mu_kk = U*b_kk;
        %gb_kk = table' * (mu_kk) - sum(exp(mu_kk));
        if(norm(d_k) <= eps) % If step to b_kk is too small
            disp(['Warning: linesearch step to small: ',num2str(norm(d_k))]);
            if(table' * (mu_kk) - sum(exp(mu_kk)) <= table' * (mu_k) - sum(exp(mu_k))) % If value of the function has not increased
                b_kk = b_k;
                disp('Exiting...');
                return;
            end 
            break;
        end
    end

    disp('Advancing:')
    ((table' * (mu_kk) - sum(exp(mu_kk)) - table' * (mu_k) + sum(exp(mu_k)) >=0))
    disp('equality?')
    disp((b_kk == b_k))
    %disp(' ')
    if(isequal(b_k,b_kk))
       disp('Convergence achieved in b_k...Exiting')
       return;
    end
    
% Update. b_kk is equivalent to  b_k + a_k .* d_k where a_k is the steepest ascending value;
    tem=b_k;
    b_k = b_kk;
    mu_k = U * b_k;
    m_k = exp(mu_k);

end

%disp('Convergence achieved:')
%disp(table' * (mu_kk) - sum(exp(mu_kk)) - table' * (mu_k) + sum(exp(mu_k)))
%disp( -table' * (mu_kk) + sum(exp(mu_kk)) + table' * (mu_k) - sum(exp(mu_k)))
%disp(b_k - b_kk)
%%disp(gb_k - gb_kk)

%%%% END 




%c_k = U*d_k;
%a_k = Get_alpha(b_k, c_k, m_k, table);

% First and second derivative at alpha_k
%phi_p_k = d'*tab - m .* exp(alpha_k .* c);
%phi_h_k = - sum( (d.^2) .* m .* exp(d .* alpha_k) );

