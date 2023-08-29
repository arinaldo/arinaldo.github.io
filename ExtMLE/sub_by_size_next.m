function [ a_new, size_new, more_new ] = sub_by_size_next ( n, a, size, more )

%% SUB_BY_SIZE_NEXT returns all subsets of an N set, in order of size.
%
%  Example:
%
%    N = 4:
%
%    1 2 3 4
%    1 2 3
%    1 2 4
%    1 3 4
%    1 3
%    1 4
%    2 3
%    1
%    2
%    3
%    (the empty set)
%
%  Discussion:
%
%    The subsets are returned in decreasing order of size, with the
%    empty set last.
%
%    For a given size K, the K subsets are returned in lexicographic order.
%
%    On the first call, it is only important that MORE be set FALSE.  The
%    input values of A and SIZE are not important.
%
%  Modified:
%
%    05 July 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the size of the set.
%
%    Input, integer A(N), the previous output subset.
%
%    Input, integer SIZE, the size of the previous output subset.
%
%    Input, logical MORE, is FALSE on the first call, which signals
%    the routine to initialize itself.  Thereafter, MORE should be TRUE.
%
%    Output, integer A_NEW(N), the next subset.
%
%    Output, integer SIZE_NEW, the size of the next subset.
%
%    Output, logical MORE_NEW, is TRUE as long as there are even more subsets
%    that can be produced by further calls.
%
  persistent more2;

  more_new = more;

  if ( ~more_new )
    more_new = 1;
    more2 = 0;
    size_new = n;
  else
    size_new = size;
    if ( ~more2 )
      size_new = size_new - 1;
    end
  end
%
%  Compute the next subset of size SIZE.
%
  if ( 0 < size_new )
    [ a_new, more2 ] = ksub_next ( n, size_new, a, more2 );
  else
    a_new = [];
    more_new = 0;
  end


  
  
  % Alternatively, to generate all the possible 0-1 patterns with N
  % elements, type:
  % >> rem(floor([0:2^N-1]' * pow2(1-N:0)),2)
  % this will produce a 2^N x N matrix.