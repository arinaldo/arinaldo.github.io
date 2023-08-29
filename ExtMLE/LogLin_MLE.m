function [MLEst, df, face_vertex, LR, Pears_chi2, LR_pvalue,Pears_chi2_pvalue] = LogLin_MLE(counts,Table,varargin)%,type,ipf) 
% LogLin_MLE: compute the Extended MLE for Log-Linear Models.
% Inputs:
%   counts: table counts.
%   Table: structure encoding the Log-Linear Model.
%   type: type of method for finding facial set. The default is Linear
%         Programming. Set it to 'N' for the Newton Method of unconstrained
%         optimization.
% Output:
%   MLEst: the Extended MLE.
%   dim: the number of parameters under Poisson sampling. It is equal to
%        the number of columns of the design matrix U used to compute the
%        Extended MLE.
%   face_vertex: set of indexes identifying the face of the convex support 
%                containing the sufficient statistics in its relative 
%                interior.


% NOTE ON ROUNDING. I use myzero

    
% Check for appropriate input arguments. 
if(nargin<2)
    error('Not enough input arguments for the function LogLin_MLE...exiting.'); 
end
nargs=length(varargin);
type_facial = 'LP';
type_ipf='N';
maxiter = 50;
epsipf = [];
maxiter = [];

if(nargs>0)
    for i=1:2:nargs        
        switch lower(varargin{i})
            case 'type'
                type_facial=varargin{i+1};            
            case 'ipf'
                type_ipf = 'Y';
            case 'maxiter'
                maxiter = varargin{i+1};
            case 'epsipf' 
                epsipf = varargin{i+1};
            case 'type_facial'
                type_facial = varargin{i+1};
            otherwise    
        end
    end
end


% Use IPF
if (type_ipf=='Y')
    if(Table.Hier=='H')
        A=Make_Basis(Table,'S');
        t = A'*counts;
        if(Is_Triangulated(Table.Complex))        
            MLEst = IPF(counts,Table,'A',A','t',t,'type','D');            
            %face_vertex = find(sum(A(:,find(t==0)),2)==0);
            face_vertex = find(MLEst>0);
            df = length(face_vertex) - rank(A(face_vertex,:));
            %dim = rank(A(face_vertex,:));
            %A = Recompute_Basis(A(face_vertex,:));
            %dim = size(A,2);
            %U = Make_Basis(Table);
            %U = Recompute_Basis(U(face_vertex,:));
            %dim = size(U,2);            
            [LR, Pears_chi2, LR_pvalue, Pears_chi2_pvalue] = Get_Chi2(counts,MLEst,df);
        else    
            % Adjust tolerance level, if necessary.            
            tol = 1e-8;
            
            MLEst=IPF(counts,Table,'A',A','t',t,'maxiter',maxiter,'eps',epsipf);
            dim = Compute_Dim(Table);
            face_vertex = 1:length(counts);
            face_vertex(MLEst<tol)=[];
            df = length(face_vertex) - rank(A(face_vertex,:));
            %dim = rank(A(face_vertex,:));
            %A = Recompute_Basis(A(face_vertex,:));
            %dim = size(A,2);
            %U = Make_Basis(Table);
            %U = Recompute_Basis(U(face_vertex,:));
            %dim = size(U,2);
            [LR, Pears_chi2, LR_pvalue, Pears_chi2_pvalue] = Get_Chi2(counts,MLEst,df);
        end
        return
    else
        error('IPF algorithm implemented only for Hierarchical Log-Linear Models...exiting.')
    end
end


% Global variable for zero.
global myzero
U = [];
face_vertex = [];
myzero=0;%1e-5; this can be set to zero since I am foring integer entries 
              % all the relevant matrices



% Check for non-existence of the MLE only if there are sampling zeroes.
if(any(counts==0))
    
    coface = [];
    first_coface=[];
    % ===================  First check for null margins. ==================
    % This first check for zero margins (which gives a facial set) is not 
    % necessary, in principle.
    % In fact, the Newton method is able to identify these cases. However,
    % it might be of some computational advantage to eliminate these 
    % easy-to-find facial sets.
    
    % Compute the 0-1 (redundant) matrix basis.
    A=Make_Basis(Table,'S');
    
    % Check for null margins.
    t=A'*counts;
    
    % If there are any, recompute the basis matrix after expunging the 
    % corresponding rows. The new U matrix is formed by a subset of the 
    % columns of the original U matrix.
    if(any(t==0))
        %face_vertex=abs(max(-1,-sum(A(:,find(t==0)),2))+1);
        face_vertex=find(sum(A(:,find(t==0)),2)==0);        
    else
        face_vertex=1:length(counts); 
    end
    % =====================================================================
    
   
       
    % ===================   Decomposable case and IPF  ====================
    % Check whether the model is decomposable, in which case use the IPF 
    % algorithm to compute the Extended MLE. 
    if(Is_Triangulated(Table.Complex)&(Table.Hier=='H'))        
        MLEst=IPF(counts,Table,'A',A','t',t,'type','D');
        df = length(face_vertex) - rank(A(face_vertex,:));
        [LR, Pears_chi2, LR_pvalue, Pears_chi2_pvalue] = Get_Chi2(counts,MLEst,df);
        return;        
    end        
    % =====================================================================
    
   
    % Compute non-sparse, full-rank basis matrix
    U=Make_Basis(Table);
    % Keep the design matrix full-rank, if necessary.
    if(length(face_vertex) < length(counts))        
        U=Recompute_Basis(U(face_vertex,:)); 
    end
    
    % ================  Check for existence of the MLE.  ==================
    % Check existence of the MLE and possiblly update face-vertex incidence
    % vector
    len_face=length(face_vertex);
    supp=find(counts(face_vertex)>0); % supp is a subset of
                                      % 1:length(face_vertex)   
        
    rank_supp = rank(U(supp,:));
    %rank_supp = -1;
    
    % First chek whether Uznot has a non-trivial null space. If it
    % doesn't, then no need to check further because face_vertex is 
    % already the face-vertex incidence vector for the face containing the 
    % observed margins in its relative interior.
        
    if((rank_supp < size(U,2)) & (length(supp) < len_face)) % U is full column rank
        % Uznot has a non-trivial null space.        
        % Second condition of the above if-statement checks that the only 
        % null cells in the table are the ones associated to zero margins.
 
        suppnot = 1:len_face;
        suppnot(supp)=[];
                   
        Uz=U(suppnot,:);
        Uznot=U(supp,:);
        % Compute an integer basis for null space of Uznot using rref.
        X = Make_Integer(null(Uznot,'r'));
        % Compute a non-integer basis for null space of Uznot using svd.
        % X=null(Uznot);

        Y = Uz * X;            

        if(size(Y,2)==1)
            % Y is a vector.

            if(all(Y>=myzero)||(all(Y<=myzero)))                                
                coface=find(abs(Y)>myzero);                                
            else
                coface=[];
            end            
        else                        
            % Y is a matrix.           
            Ysum=sum(Y,2);    
            Yplus=[];
            Yzero=[];
                
            % Before looking for facial sets, do simple checks.
            % First easy check. 
            if(all(Ysum==myzero))                  
                % The MLE exists by Stiemke's Theorem
                % of Alternatives.
                coface=[]; 
            else                                             
                % Checks
                [Yzero,Yplus]= Preprocess(Y,Ysum);
                if(length(Yplus)>0)
                   first_coface =  suppnot(Yplus);
                end
                    
                % Decide which rows to keep. They correspond to the union 
                % of the index vectors Ypositivesum,Ysamesign and Yzero 
                % (which is disjoint from teh others).                    
                %length(Ypositivesum)
                if(length(Yplus) + length(Yzero) == length(suppnot)) 
                    % The co-face is found.
                    coface=first_coface;
                else                     
                    % The co-face must be found with other methods.
                    Ykeeprows = 1:size(Y,1);
                    Ykeeprows([Yplus;Yzero]) = [];
                    Y=Y(Ykeeprows,:);
                    suppnot=suppnot(Ykeeprows);
                    
                    % If Y not of full-column rank, remove redundant 
                    % columns. 
                    if(rank(Y)<size(Y,2))                        
                        Y = Recompute_Basis(Y);                          
                    end

                    fprintf('\n===============================================================================\nComputing the facial set using')
                    
                    % Compute the vector coface. It is a subset of
                    switch type_facial
                        case 'LP'                            
                            fprintf(' LP methods...\n\n')
                            [coface, exitflag]=Find_Facial_Set(Y,type_facial);
                            switch exitflag
                                case 0
                                    [coface, exitflag]=Compute_Facial_LP(Y,type_facial,300);
                                    switch exitflag
                                        case 0
                                             error('Error in Compute_Facial_LP (MaxIter). Exiting'); 
                                        case {-2,-3,-4,-5}
                                            error('Error in Compute_Facial_LP. Exiting'); 
                                    end                                                                                                
                                case {-2,-3,-4,-5}
                                    error('Error in Compute_Facial_LP. Exiting');   
                                otherwise                                
                            end
                        case 'N'
                            fprintf(' Newton method...\n\n')    
                            [coface,exitflag]=Find_Facial_Set(Y,type_facial);
                            switch exitflag
                                case 0
                                    [coface,exitflag]=Find_Facial_Set(Y,type_facial,300);
                                    switch exitflag
                                        case 0
                                             error('Error in Compute_Facial_Newton (MaxIter). Exiting'); 
                                        case {-1}
                                            error('Error in Compute_Facial_Newton. Exiting'); 
                                    end                                                                                                             
                                case -1
                                    error('Error in Compute_Facial_Newton. Exiting');                                                                
                                otherwise                                
                            end                           
                        case 'E1'
                            fprintf(' Entropy methods.\n\n');
                            tol = 1e-8;
                            [coface, exitflag] = Find_Facial_Set(Y,type_facial,tol);                            
                        case 'E2'
                            fprintf(' Entropy methods.\n\n');
                            tol = 1e-8;
                            [coface, exitflag] = Find_Facial_Set(Y,type_facial,tol);                                                        
                        otherwise                            
                            fprintf(' Newton method (Illegal option %s. Using Newton method instead...\n\n',type_facial);
                            [coface, exitflag] = Compute_Facial_Newton(Y);                        
                    end
                    
                    fprintf('\n...done!\n')
                    
                end    
            end
        end                                   
    end   
    
    
    
    % If coface and/or first _coface is not empty, recompute the basis 
    % (design) matrix as above.
    if( (length(coface)>0)||(length(first_coface)>0) )
        face_vertex = face_vertex( setdiff( (1:len_face), [suppnot(coface);first_coface]) );
        U=Recompute_Basis( U(setdiff( (1:len_face), [suppnot(coface);first_coface] ),:) );         
    end    
    %======================================================================    
end 

if(isempty(U))
   U = Make_Basis(Table); 
end
if(isempty(face_vertex))
   face_vertex = 1:size(U,1); 
end

fprintf('\n===============================================================================\nComputing the MLE...\n\n');
theta_hat=Newton(U,counts(face_vertex));
%theta_hat=Matlab_Newton(U,counts(face_vertex));
fprintf('\n...done!\n\n');

MLEst=zeros(length(counts), 1);
MLEst(face_vertex)=exp(U*theta_hat);

% Compute the appropriate dimension of the parameter space. It is possible
% to use Cholesky factorization instead, which is to be done in a 
% non-Matlab implementation.
df = length(face_vertex) - size(U,2);
[LR, Pears_chi2, LR_pvalue, Pears_chi2_pvalue] = Get_Chi2(counts,MLEst,df);





%%%% SUBFUNCTIONS


function [coface,exitflag] = Find_Facial_Set(B,routine,maxiter)

if(nargin==2)
    maxiter=[];
end

switch routine
    case 'N'
        % Newton method
        if(length(maxiter)==0)
            [coface,exitflag]=Compute_Facial_Newton(B);
        else
            [coface,exitflag]=Compute_Facial_Newton(B,maxiter);
        end        
    case 'LP'        
       % LP method
        if(length(maxiter)==0)
            [coface, exitflag]=Compute_Facial_LP(B);
        else
             [coface, exitflag]=Compute_Facial_LP(B,maxiter);
        end
    case 'E1'        
        % Entropy 1
        tol = 1e-8;
        [coface, exitflag] = Compute_Facial_Entropy1(B,tol);
    case 'E2'
        % Entropy 2
        tol = 1e-8;        
        [coface, exitflag] = Compute_Facial_Entropy2(B,tol);
    otherwise        
end







function [lr, pears_chi2, lr_pvalue, pears_chi2_pvalue] = Get_Chi2(x,xmle,df)
% Compute Likelihood ratio, Pearson chi^2 and their pvalues
lr = 2 * sum( x(x>0) .* log( x(x>0)./xmle(x>0) ) );
pears_chi2 = sum ( ( (x(xmle>0) - xmle(xmle>0) ).^2 ) ./ xmle(xmle>0) );
pvalues = 1 - chi2cdf([lr pears_chi2],df);
lr_pvalue = pvalues(1);
pears_chi2_pvalue = pvalues(2); 





function [Xzero,Xplus]= Preprocess(X,Xsum)
% Compute the index vectors: Xpositivesum, Xsamesign and Xzero.
Xpositivesum=[];
if(nargin < 2)
   Xsum = sum(X,2); 
end
myzero = 0;

if( (all(Xsum >= myzero))||(all(Xsum <= -myzero)) ) 
% Check whether the rows of X sum to a positive or 
% negative vector. If so, the non-zero coordinates are 
% part of the co-face.
% If all the coordinatesof Xsum are non-zero and of the
% same sign, then the zero cell define a co-face.
Xpositivesum=find(sum(abs(X),2)>myzero);                         
end

% Check for columns having the same sign and eliminate the corresponding 
% rows.    
Xsamesign=find(all(X>=myzero,1)+all(X<=-myzero,1)>0);                    
if(length(Xsamesign)>0)    
    Xsamesign = find(sum(abs(X(:,Xsamesign)),2)>myzero);                          
end
    
% Compute the appropriate indices.

% Check whether there are null rows of Y.
Xzero=find(sum(abs(X),2)==myzero);
Xplus = sort(unique( [Xpositivesum; Xsamesign] ));

% Repeat check using rref     
if(length(Xplus) + length(Xzero) < size(X,1))      

    XXpositivesum=[];
    keep = (1:size(X,1))';
    keep([Xplus;Xzero]) = [];
    XX =  Make_Integer(rref(X(keep,:)')');
    XXsum = sum(XX,2);     
    if( (all(XXsum >= myzero))||(all(XXsum <= -myzero)) ) 
        XXpositivesum=find(sum(abs(XX),2)>myzero);                         
    end
    XXsamesign=find(all(XX>=myzero,1)+all(XX<=-myzero,1)>0);                    
    if(length(XXsamesign)>0)    
        XXsamesign = find(sum(abs(XX(:,XXsamesign)),2)>myzero);                          
    end
    XXplus = sort(unique( [XXpositivesum; XXsamesign] ));
    Xplus = sort(unique([Xplus;keep(XXplus)]));
    
end
                    
                
    

                   