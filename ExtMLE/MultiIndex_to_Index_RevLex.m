function index = MultiIndex_to_Index_RevLex(siz,multiindex)
% Same as sub2ind, except it accepts a vector of coordinates as a second
% argument.
% Inputs:
%   multiindex: multiindex value.
%   siz: size vector.
% Outpus:
%   multiindex: multiindex

if (length(siz)<2)
        error('Size vector must have at least 2 elements.');
end

if( (length(multiindex)~=length(siz)) || (prod(multiindex)> prod(siz)))
   error('Input vectors must be compatible');
elseif (sum(multiindex > siz)>0)   
    error('First input vector >= second input vector.');
end
    

prodind = [1 cumprod(siz(1:end-1))];
index = 1 + sum((multiindex-1).*prodind);
