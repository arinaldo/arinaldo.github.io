function multiindex = Index_to_MultiIndex_RevLex(index,siz)
% Index to MultiIndex using revlex order
% Inputs:
%   index: index value.
%   siz: size vector
% Outputs:
%   multiindex: multiindex.

if (length(siz)<2)
        error('Size vector must have at least 2 elements.');
end

if(prod(siz)<index)
   error('Index exceed size vector. Exiting...'); 
end

prodsize = cumprod(siz(1:(end-1)));
lensize = length(siz);
multiindex=zeros(lensize,1);

for i=lensize:-1:1
    if (i==1)
        multiindex(i) = index;
    else
        multiindex(i) = floor((index-1)/prodsize(i-1)) + 1;        
        index = index - prodsize(i-1)*(multiindex(i)-1);
    end
end