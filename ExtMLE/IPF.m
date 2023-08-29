function m = IPF(counts,Table,varargin)%IPF(x,table,type,A,eps,maxiter)
% IPF: Iterative Proportional Fitting algorithm.
%   See "Graphical Models" by Lauritzen.
% Inputs:
%   counts: table counts.
%   Table: structure encoding the Log-Linear Model.
%   varargin (optional arguments): 
%       A: sparse design matrix. Can be computed from the structure Table.
%       t: margins, computed as A*counts.
%       type: type of Log-Linear Model. If equal to 'D' it specifies a
%             decomposable model. Default: not decomposable.
%       eps: maximum tolerance between fiited and observed margins. 
%            Default = 1e-8.
%       maxiter: maximum number of IPF cycles. Default = 50.
% Ouputs:
%   m: Maximum Likelihood Estimates.


% =========== Input arguments and defaulting of the parameters. ===========
if (nargin<2)
   error('Error in function IPF: not enough arguments. Abort.') 
end


eps=1e-8;
maxiter=50;
type='ND';

nargs=length(varargin);
if(nargs>0)
    % Sparse design matrix.
    Ai=[];
    for i =1:2:nargs
        if(varargin{i}=='A')            
            Ai = i;
            A=varargin{Ai+1};           
            break
        end        
    end
    if(isempty(Ai)), A=Make_Basis(Table,'S')'; end
    
    % Margins.
    ti=[];
    for i =1:2:nargs
        if(varargin{i}=='t')            
            ti = i;
            t=varargin{ti+1};           
            break
        end        
    end
    if(isempty(ti)), t=A*counts; end
    
    % Other arguments.
    for i =1:2:nargs               
        switch varargin{i}        
            case 'type'                        
                type=varargin{i+1};
            case 'A'  
            case 't'
            case 'eps'
                if(~isempty(varargin{i+1})), eps=varargin{i+1}; end;
            case 'maxiter'
                if(~isempty(varargin{i+1})), maxiter = varargin{i+1}; end;
            otherwise
                disp(['WARNING: ' varargin{i} ' is not a valid argument for the function IPF. ']);
        end
    end
else
    A=Make_Basis(Table,'S')';
    t=A*counts;
end
% =========================================================================



% ========================   IPF Algorithm.   =============================
% Initialize.
m=ones(length(counts),1);
mt=A*m;    
iter=1;
 if(type=='D')
        % If the model is decomposable, then it is possible to obtain
        % convergence in just one iteration by using a perfect ordering
        % return; 
    end
while((iter<maxiter)&(max(abs(t-mt)) > eps))
   % IPF cycle.
    k=0;
    for a=1:length(Table.Complex)
        nia=prod(Table.Ncat(Table.Complex{a}));
        Aa=A((k+1):(k+nia),:);
        ma = Aa*m;
        ta=t((k+1):(k+nia));
        adjust=zeros(length(ta),1);
        adjust(ma>0) = ta(ma>0) ./ ma(ma>0);
        Aa_adjust = diag(adjust) * Aa;
        m = m .* sum(Aa_adjust,1)';     
        k=k+nia;
    end 
    iter=iter+1;
    mt = A*m;
end
% =========================================================================

if(max(abs(t-mt)) > eps)               
    disp('================= IPF WARNING ==================');       
    disp('========== Convergence not achieved. ===========');
    disp(['Maximum number of iterations: ' num2str(maxiter) '.']);
    disp(['Default tolerance: ' num2str(eps) '.']);
    disp(['Achieved deviance: ' num2str(max(abs(t-mt))) '.']);
    disp('================================================');
end
    
