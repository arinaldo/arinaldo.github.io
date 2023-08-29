function U = Make_Basis(Table,type,complex)
% Make_Basis: build the design matrix.
%   Number of rows is the numebr of cells and number of columns the number
%   of parameters under Poisson sampling scheme.
% Inputs: 
%   Table: structure encoding the Log-Linear Model.
%   type: type of design matrix. The default = 'C', for the contrast basis. 
%         If type = 'S', then the sparse 0-1 design matrix is computed.
%   complex: class of subsets holding the interactions. 
% Inputs:
%   U: the design matrix. size(U,1) = number of cells,
%                         size(u,2) = number of parameters under Poisson
%                                     sampling.

if(nargin<2)
    type='C';
end

% Default for complex
if (nargin<3)
    complex=Table.Complex;
end

if(type=='C')
    if(Table.Hier=='H')
       Faces=Compute_Faces(Table.Complex, Table.Nvar); 
    else
        Faces=complex;
    end
    for i=1:length(Faces)        
        if(Faces{i}==0) 
           Ui=ones(prod(Table.Ncat),1); 
        else
            Ui=1;
            for j=1:Table.Nvar
                if(ismember(j,Faces{i}))
                   Z=eye(Table.Ncat(j)) +  diag(-ones([Table.Ncat(j)-1, 1]),-1);
                   Z=Z(:,(1:Table.Ncat(j)-1));
                   % Different way of generating Z:
                   %Z = [ones(1,Table.Ncat(j)-1); -eye(Table.Ncat(j)-1)];
                else
                    Z=ones([Table.Ncat(j) 1]);
                end;
                Ui=kron(Ui,Z);                
            end
        end
        if(i==1)
            U=Ui;
        else
            U=[U Ui];
        end
    end
else
    if(Table.Hier ~= 'H')
        error('ERROR: Sparse basis is defined only for a Hierarchical Model. Exiting.');
    end
    Facets=Table.Complex;
    for i=1:length(Facets)        
        Wi=1;
        for j=1:Table.Nvar            
            if(ismember(j,Facets{i})) 
               W=eye(Table.Ncat(j));
            else
               W=ones([Table.Ncat(j) 1]);
            end
            Wi=kron(Wi,W);
        end
        if(i==1)
            U=Wi;
        else
            U=[U Wi];
        end        
    end
end
   







