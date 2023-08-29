function graphical = Is_Graphical(Complex, G)
% Is_Graphical: check whether a Log-Linear model encoded by the simplicial
%   complex Complex is graphical or not. 
%   The clique hypergraph of the interaction graphs contains a hyper-edge 
%   which is strictly bigger than some facet of the simplicial complex 
%   if and only if the model is not graphical.
%   See "Graphical Models" by Lauritzen, page 89.
% Inputs:
%   Complex: simplicial complex/hypergraph encoding the Log-Linear Model.
%   G:       vertex-incidence matrix for the associated interaction graph 
%            (optional).
% Outputs:
%   graphical: 1 if the model is graphical and 0 otherwise.


graphical=1;
ncliques = length(Complex);

% If G matrix is not provided.
if(nargin<2)    
    % Compute the vertex-incidence matrix for the interaction graph.    
    % Allocate memory for G
    dimG = max(Complex{1});
    for i = 2:ncliques
        dimG = max(max(Complex{i}),dimG);
    end
    G = zeros(dimG);
    %
    for j=1:ncliques
        G(Complex{j},Complex{j}) = 1; 
    end
end

nvertices = length(G);

% Check whether the model is graphical.
for j=1:ncliques
    noncliquej=setdiff( 1:nvertices, Complex{j} );
    for i=1:length(noncliquej)
      Gindex=sort([noncliquej(i);Complex{j}]);
      lGindex=length(Gindex);
      if( sum(G(Gindex,Gindex),2) == lGindex .* ones(lGindex,1) )
          graphical=0;
          return;
      end
    end
end





%%=========================================================================
%% This version of the function Is_Graphical(Complex, G)will actually find 
%% the whole interaction graph. 
%% WARNING: NOT TESTED!

%[Gu,Iu]=unique((G+eye(nvertices))','rows');
%k=1;
%Iu=sort(Iu);
%for g=1:length(Iu)    
%    bdg=find(G(:,Iu(g))==1); % boundary of vertex Iu(g), sorted.
%    % check for complete subgraph induced by bdg
%    nbdg=length(bdg);
%    %subGcomplete=[]; % complete subgraph of G
%    subGtochek=dbg; % sugraph of G to be checked for completeness
%    
%    while(length(subGtochek)>0)
%        [subGcomplete,subGtochek]=Check_Complete_SubGrpah(G,subGtochek);
%        temp=sort(subGcomplete);
%        if(k==1)            
%            new_complex{k}=temp;
%            k=k+1;            
%        else
%            is_new=1;
%            % Check whether the proposed complete subgraph is not already
%            % present in the simplicial complex/hpergraph
%            for j=1:(k-1)
%                if( (length(new_complex{j}) >= length(temp)) & (sum(ismember(temp,new_comple{j}))==length(temp)) ) %% check this               
%                    is_new=0;
%                    break;
%                end
%            end
%       end
%       if(is_new)           
%           new_complex{k}=temp;
%           k=k+1;
%       end     
%       
%    end %  while(length(subGtochek)>0)     
%   
%end %for g=1:length(Iu)    




%function [complete,tocheck]=Check_Complete_SubGrpah(Graph,subtocheck)
% 
%a_new=[];
%size_new=0;
%more_new=0;
%N=length(subtocheck);
%while((size_new>2)||(size_new==0))
%    
%    [a_new, size_new, more_new] = sub_by_size_next(N,a_new,size_new,more_new);
%    % Check for a complete subgraph.
%    if(size_new>2)
%        index=subtocheck(a_new);
%        if(isequal(Graph(index,index),ones(size_new)))
%           complete =  subtocheck(a_new);
%           tocheck = setdiff(subtocheck,complete);
%           return;
%        end
%    else%if(size_new==2)        
%        complete =  subtocheck(a_new);
%        tocheck = setdiff(subtocheck,complete);
%        return;
%    end    
%end
%%=========================================================================







