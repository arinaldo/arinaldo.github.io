function multiindex = Index_to_MultiIndex_Lex(index,siz)
% Index to MultiIndex using lex order
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


prodsize = fliplr(cumprod(siz(end:-1:2)));
lensize = length(siz);
multiindex=zeros(lensize,1);

for i=1:lensize
    if (i==lensize)
        multiindex(i) = index;
    else
        multiindex(i) = floor((index-1)/prodsize(i))+1;  
        index = index - prodsize(i)*(multiindex(i)-1);
    end
end