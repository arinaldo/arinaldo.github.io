function [triangulated, labels] = Is_Triangulated(Complex)
% Is_Triangulated: perform the Maximum Cardinality Search Algorithm.
%   The inputs must correspond to a graphical model. It uses the function
%   Is_Traingulated to check if this is the case.
%   See "Probabilistic Network and Expert Systems" by Cowell, Dawid,
%   Lauritzen and Spiegelhalter, pages 55-56.
%   WARNING: the cells of Complex are column vectors.
%
% Inputs:
%   Complex: a simplicial complex or hypergraph containing the cliques of 
%            the underlying graph. It is a cell array. It is assumed that 
%            the vertexs are labelled from 1 to n, the number of vertices.
% Output:
%  triangulated: 1 if the underlying graph is triangulated and 0 otherwise.  
%  labels: a perfect numbering if the graph is triangulated or an empty 
%          vector otherwise.


% Compute the vertex-incidence matrix for the interaction graph.
ncliques = length(Complex);
for j=1:ncliques
   G(Complex{j},Complex{j}) = 1; 
end
nvertices = length(G);


% First Check whether Complex is associated to a graphical model.
if(~Is_Graphical(Complex,G))
   triangulated=0;
   labels=[];
   return;
end

% Proceed with the Maximum Cardinality Search Algorithm.
triangulated=1;
labels=[];
i=1;

% Vertex index vector and weight vector.
V=1:nvertices;
weights=zeros(nvertices,1);

L=[];

while(~isequal(L,V))    
    U=setdiff(V,L);
    if(isempty(U)), return; end;
    [maxv,maxvindex]=max(weights(U));
    v=U(maxvindex);
    labels(v)=i;
    ne = Neighbors(G(:,v),v);
    neL = intersect(ne,L);
    
    if( ~isequal(G(neL,neL), ones(length(neL)) ))
       triangulated = 0;
       labels=[];
       return;
    end
    
    neU=intersect(ne,U);
    weights(neU) =  weights(neU) + 1;
    L=[L v];
    i=i+1;
end



function nb = Neighbors(Gv,v)
% Function that returns the neighbors of the vertex v in the graph G. G is
% encoded as a 0-1 vertex-incidence symmatrix matrix. 
% The output is a column vector.
Gv(v)=0;
% This vector contains the neighbors of v.
nb = find(Gv==1);











