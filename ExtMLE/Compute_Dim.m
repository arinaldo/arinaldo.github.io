function dim = Compute_Dim(Table, varargin)%type, Faces)
% Compute_Dim: compute the dimension of the Log-Linear model specified by
%   the cell-array Table.
% Inputs:
%   Table: structure encoding the Log-Linear Model.
%   type: type of Log-Linear Model. default = 'H', for a Hierarchical
%         Model.
%   Facets: Cell arrays representing a hypergrgaph on Table.Nvar nodes. 
% Outputs:
%   dim: dimension of the corresponding Log-Linear Model under Poisson
%        sampling scheme.



% Set the parameters.
nargs=length(varargin);
if(nargs>0)
    for i=1:2:nargs        
        switch lower(varargin{i})
            case 'type'
                type_model=varargin{i+1};    
            case 'facets'
                facets = varargin{i+1};
        end
    end
else
    type = Table.Hier;
    facets = Table.Complex;
end

if(type=='H')
    Faces = Compute_Faces(facets, Table.Nvar);
else
    Faces = facets;
end


% Compute the dimension.
nfaces=length(Faces);
dim=0;
for i=1:nfaces
    if(Faces{i}==0)
        dim=1;
    else
        dim=dim+prod(Table.Ncat(Faces{i})-1);        
    end
end