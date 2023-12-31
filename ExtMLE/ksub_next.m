function [ a_new, more_new ] = ksub_next ( n, k, a, more )

%% KSUB_NEXT generates the K subsets of an N set, one at a time.
%
%  Modified:
%
%    06 August 2004
%
%  Reference:
%
%    A Nijenhuis and H Wilf,
%    Combinatorial Algorithms,
%    Academic Press, 1978, second edition,
%    ISBN 0-12-519260-6.
%
%  Parameters:
%
%    Input, integer N, the size of the set from which subsets are drawn.
%
%    Input, integer K, the desired size of the subsets.  K must
%    be between 0 and N.
%
%    Input, integer A(K).  A(I) is the I-th element of the
%    subset.  Thus A(I) will be an integer between 1 and N.
%    Note that the routine will return the values in A
%    in sorted order: 1 <= A(1) < A(2) < ... < A(K) <= N
%
%    Input, logical MORE.  Set MORE = FALSE before first call
%    for a new sequence of subsets.  It then is set and remains
%    TRUE as long as the subset computed on this call is not the
%    final one.  When the final subset is computed, MORE is set to
%    FALSE as a signal that the computation is done.
%
%    Output, integer A_NEW(K).  A(I) is the I-th element of the
%    subset.  Thus A(I) will be an integer between 1 and N.
%    Note that the routine will return the values in A
%    in sorted order: 1 <= A(1) < A(2) < ... < A(K) <= N
%
%    Output, logical MORE_NEW.  Set MORE = FALSE before first call
%    for a new sequence of subsets.  It then is set and remains
%    TRUE as long as the subset computed on this call is not the
%    final one.  When the final subset is computed, MORE is set to
%    FALSE as a signal that the computation is done.
%
  persistent m;
  persistent m2;

  if ( k < 0 | n < k )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'KSUB_NEXT - Fatal error!\n' );
    fprintf ( 1, '  N = %d\n', n );
    fprintf ( 1, '  K = %d\n', k );
    fprintf ( 1, '  but 0 <= K <= N is required!\n' );
    error ( 'KSUB_NEXT - Fatal error!' );
  end

  more_new = more;

  if ( ~more )
    m2 = 0;
    m = k;
  else
    a_new(1:k) = a(1:k);
    if ( m2 < n-m )
      m = 0;
    end
    m = m + 1;
    m2 = a_new(k+1-m);
  end

  for ( j = 1 : m )
    a_new(k+j-m) = m2 + j;
  end

  more_new = ( a_new(1) ~= (n-k+1) );
