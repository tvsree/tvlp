% ----------------------------------------------------------------
% Evaluate a Chebyshev series through recursion; (for LPC2LSF)
% 
% To call: y = chebyshev(gp,n,x); gp(1) to gp(n) are the
%              coefficients of the series corresponding to 
%              Chebyshev polynomials of order n-1 to 0.
% ----------------------------------------------------------------

  function [y] = chebyshev(gp,n,x);
  
  T0 = 1; T1 = x;		% coeffs evaluated in ascending order
  y = gp(n-1)*T1 + gp(n)*T0;	% from the last, backwards.
  for j = n-2:-1:1,
      T = 2*x*T1 - T0;
      y = y + T*gp(j);
      T0 = T1;
      T1 = T;
  end
  
% .................................................................