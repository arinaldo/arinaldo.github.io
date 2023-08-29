function [coface, exitflag, output] = Compute_Facial_LP(B,naxiter)
% Compute_Facial_LP: compute the facial set using LP.
%   B is a integer matrix such that the MLE exists if and only if the
%   system Bx \gneq 0 has a solution. The co-face corresponds to a
%   solution with maximal support.
%   This linear program solves the problem:
%
%   min 1'Bx
%   s.t. Bx <= 1
%        Bx >= 0
%   
%   where B = Uz * X, with X being a matrix whose columns (integer vectors) 
%   form a basis for kernel(Uznot). Zero rows from B are already expunged.
%   The LP procedure will in general not produce the final solution if
%   applied just once. It will need to be repeated until the objective
%   value is zero or the zero cells precisely defines a coface.
%   Rounding to zero of small elements of Bx is performed in 
%   order to compute the facial set.
% Inputs:
%   B: matrix such that size(B,2) = number of parameters under Poisson
%                                   sampling.
% Outputs:
%   coface: index set containing the coface with maximal support.
%           coface is a subset of 1:size(B,1).
%
% Comments: the following options are used for linprog:
%           optimset('Simplex','off','LargeScale','on','Display','iter')
% so that interior points method is used instead of the simplex method.

fval=-1;
firsttime=1;
myzero=1e-5;%-1e-16;
if(nargin==1)
    options=optimset('Simplex','off','LargeScale','on','Display','iter');
else
    options=optimset('Simplex','off','LargeScale','on','Display','iter','MaxIter',maxiter);
end
ncycle=1;
allrB = size(B,1);


% Cycles.
while(fval < -myzero) 
    if(firsttime)
        fprintf('Performing cycle %5d ...', ncycle);       
        % Setting up the variables for linprog
        rB=allrB;
        cB=size(B,2);        
        f=-ones(1,rB)*B;        
        A=[B;-B];
        b=[ones(rB,1).*10;zeros(rB,1)];

        % Calling the LP minimization routine.
        [x,fval,exitflag,output]=linprog(f,A,b,[],[],[],[],[],options);            
        %=====================================
        y = B*x;
        if(abs(fval)<=myzero)
            coface=[];
            return;
        else                  
            next_coface=find(y<= myzero); 
            if(length(next_coface)<=1)
                coface=1:allrB;
                return;
            end
        end
        firsttime=0;
        ncycle=ncycle+1;
    else        
        fprintf('Performing cycle %06d ...\n', ncycle);
        disp(strcat('Performing cycle ',num2str(ncycle),'...'));
        % Re-computing matrix B and linprog variables after expunging coface.
        B= Recompute_Basis(B(next_coface,:));               
        rB=size(B,1);
        cB=size(B,2);
        f=-ones(1,rB)*B;        
        A=[B;-B];
        b=[ones(rB,1).*10;zeros(rB,1)];           

        % Calling the LP minimization routine.
        [x,fval,exitflag,output]=linprog(f,A,b,[],[],[],[],[],options);
        %=====================================
        y = B*x;
       % if (y >= 1 - myzero)
        if(abs(fval)<=myzero) 
            coface = 1:allrB;
            coface(next_coface)=[];
            return;
        else                       
            next_coface=next_coface(find(y<= myzero));  
            if(length(next_coface)<=1)
                coface=1:allrB;
                return;
            end            
        end            
        ncycle=ncycle+1;                        
    end
end

