function hyper = Simulate_Hypergraph(N,p)
% Simulate_Hypergraph: simulate hypergraphs on N nodes and hence Log-Linear
%   Models.
% Inputs:
%   N: number of nodes.
%   p: if all nodes are covered but there can be more hyper-edges, continue 
%      generating hyperedges with probability p and stop with probability 
%      (1-p)
% Ouputs:
%   hyper: a hypergraph in the form of a cell array of sorted subsets of
%          1:N.

if (nargin<2)
	p=0.5;
end
subsets = combn([0 1],N);
subsets(1,:)=[];
nodes = 1:N;
hyper={};
k=1;
nodes_covered=[];
allnodes=0;

rand('state',sum(100*clock));

while(true)
    i = unidrnd(size(subsets,1));
   	hyper{k} = nodes(logical(subsets(i,:)))';
    nodes_covered = sort(unique([nodes_covered hyper{k}']));
    k=k+1;
    if(~allnodes)
        if(isequal(nodes_covered, nodes)), allnodes=1; end;
    end
    
    if(allnodes)
        num = rand;
        if(rand>=p)       
            break;            
        end 
    end
    
    subsets= Eliminate(subsets,i);  
    if (size(subsets,1)==0)
        %disp('zero')
        break        
    end;    
end


function out = Eliminate(A,i)
ai = A(i,:);
index_to_eliminate = [find(sum(A(:,~logical(ai)),2)==0);find(sum(A(:,logical(ai)),2) == sum(ai))];
out = A(setdiff(1:size(A,1),index_to_eliminate),:);

