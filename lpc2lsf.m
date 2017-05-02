% -------------------------------------------------------------
% Conversion of LPC coefficients to LSF representation using
% Chebyshev polynomial recursion. (Ref: Kabal & Ramachandran,
% IEEE Tran. ASSP, Vo.34, No.6, Dec. 1986)
% 
% To call: lsf = lpc2lsf(lpc,p); lpc(1)-lpc(p+1) are the LPCs
%          and lsf(1)-lsf(p) are the LSFs in the range of 0-pi.
% -------------------------------------------------------------
  function [lsf] = lpc2lsf(lpc,p);
% --------------------------------
  
  if (lpc(1) ~= 1),
      disp('lpc2lsf: Ak(1) ~= 1');
      return;
  end
  A = lpc(1:p+1);
  fa = [A 0]+[0 fliplr(A)];		% p+1 order aux. polynomials
  fb = [A 0]-[0 fliplr(A)];
  
  even = 0;
  if (fix(p/2) == p/2), even = 1;  end
  if even,
      na = p/2 + 1;
      nb = na;
  else
      nb = (p+1)/2;
      na = nb+1;
  end
  
  if even,
      ga(1) = fa(1);
      gb(1) = fb(1);
      for j = 2:na,		% 0-(na-1) is the order of polynomial
	  ga(j) = fa(j)-ga(j-1);
	  gb(j) = fb(j)+gb(j-1);
      end
  else
      ga = fa;
      gb = fb;
      for j = 3:nb,		% 0-(nb-1) is the order of polynomial
	  gb(j) = fb(j)+gb(j-2);
      end
  end
  for j = 1:na-1,		% chebyshev series coeffs: DECREASING order
      gpa(j) = 2*ga(j);
  end
  for j = 1:nb-1,
      gpb(j) = 2*gb(j);
  end
  gpa(na) = ga(na);
  gpb(nb) = gb(nb);

% Evaluate the series using recursive defn of Tm(x) = 2x.Tm-1(x)-Tm-2(x)
  
  delta = 0.02;			% coarse grid on x: delta < min[x(i)-x(i-2)];
  eps = 0.0015;			% root accuracy: eps < min[x(i)-x(i-1)];

  found = 0;
  nlsf = 0;
  polya = 1;			% flag indicates polynomial evaluated
  xl = 1;
  yl = chebyshev(gpa,na,1);
  x = xl-delta;
  while (x > -1 & nlsf < p),
      if (polya), 
          y = chebyshev(gpa,na,x);
      else
	  y = chebyshev(gpb,nb,x);
      end
      if (y*yl < 0),
	  for j = 1:4,
	      zcr = (xl-x)*y/(y-yl);
	      if (zcr < eps/2),		% root is FOUND!
		  found = 1;
		  x0 = x+zcr;
	      elseif (xl-x-zcr < eps/2),	% root is FOUND!
                  found = 1;
		  x0 = xl;
	      else
	          if (polya), 
		      ynew = chebyshev(gpa,na,x+zcr);
	          else
		      ynew = chebyshev(gpb,nb,x+zcr);
	          end
	          if (y*ynew < 0),
		      yl = ynew;
		      xl = x+zcr;
	          else
		      y = ynew;
		      x = x+zcr;
	          end
	          if (xl-x < eps),
		      x0 = (x+xl)/2;	% root is FOUND!
		  end
	      end 
              if (found),
		  found = 0;
	          nlsf = nlsf+1;
		  lsf(nlsf) = acos(x0);
		  polya = 1-polya;	% switch polynomials
		  x = x0;
		  if (polya),
		      y = chebyshev(gpa,na,x);
		  else
		      y = chebyshev(gpb,nb,x);
		  end
		  break;		% exit root refinement loop
	      end
	  end
      end
      yl = y;
      xl = x;
      x = x-delta;
  end
  if (nlsf < p), fprintf('lpc2lsf: ??Less roots: %3.0f\n',nlsf); end

% ........................................................................