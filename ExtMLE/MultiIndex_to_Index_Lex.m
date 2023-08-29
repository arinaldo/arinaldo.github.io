function index = MultiIndex_to_Index_Lex(siz,multiindex)
% Same as sub2ind, except it accepts a vector of coordinates as a second
% argument AND order the array entries in lexicographi order rather than
% reverse lexicographic order, which is the default order in Matlab. 
% Index to MultiIndex using lex order
% Inputs:
%   multiindex: multiindex value.
%   siz: size vector.
% Outpus:
%   index: index
if (length(siz)<2)
        error('Size vector must have at least 2 elements.');
end

if( (length(multiindex)~=length(siz)) || (prod(multiindex)> prod(siz)))
   error('Input vectors must be compatible');
elseif (sum(multiindex > siz)>0)   
    error('First input vector >= second input vector.');
end

prodind = cumprod(siz(2:end));
prodind=[prodind(end:-1:1) 1];
index = 1 + sum((multiindex-1).*prodind);
