% -----------------------------------------------------------
% Conversion of LSF parameters to LPC parameters by direct
% multiplication of 2nd order sections for each LSF.
% 
% To call: lpc = lsf2lpc(lsf,p); lsf(1)-lsf(p) are the LSFs in
%          the range of 0 to pi; the output is a stable FIR 
%          filter of lpc(1) to lpc(p+1), with lpc(1)=1;
% ------------------------------------------------------------

  function [lpc] = lsf2lpc(lsf,p);
  
% --------------------------------
  
  if (p < 1), return; end
  na = 1;				% initialization
  nb = 1;
  ga = 1;
  gb = 1;
  polya = 1;
  for j = 1:p,
      b1 = -2*cos(lsf(j));		% second order section with b2=1
      if (polya),
	  ga = [ga 0 0] + [0 b1*ga 0] + [0 0 ga];
	  na = na+2;
	  polya = 1-polya;
      else
	  gb = [gb 0 0] + [0 b1*gb 0] + [0 0 gb];
	  nb = nb+2;
	  polya = 1-polya;
      end
  end

  if (fix(p/2) == p/2),			% reconstruct aux. polynomials
      fa = [ga 0]+[0 ga];
      fb = [gb 0]-[0 gb];
  else
      fa = ga;
      fb = [gb 0 0]-[0 0 gb];
  end
  
  lpc = (fa+fb)/2;
  lpc = lpc(1:p+1);
  
% ---------------------------------------------------------------------