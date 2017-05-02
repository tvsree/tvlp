%------------------------------------------------------------
% Time-varying linear prediction coefficient (TVLPC) analysis
% of speech segments; TVLPC solution formulation
% is based on polynomial or other basis functions of LPCs. 
% The polynomial/bases order and the predictor order
% can be chosen based on segment duration and signal BW
% (sampling frequency), respectively.
%
% The least squares formulation is similar to the covariance
% formulaton of linear prediction, hence no explicit window
% function is used. 
%------------------------------------------------------------
  function [tvLPC tvLPres] = TVLPC(sig,lpcRdr,polRdr,type);
%----------------------------------------------------------
          % construct basis functions

          randPhase = [0 rand(1,polRdr)*2*pi];      
          for i = 0:polRdr,     % polynomial bases upto polRdr
              if ~isempty(strfind(type,'pol')), 
                  bases(i+1,:) = ([0:sigL-1]/sigL).^i; % polynomial bases
              end
              if ~isempty(strfind(type,'sin')), 
                  bases(i+1,:) = cos((i*2*pi*[0:sigL-1]/sigL) + randPhase(i+1)); % sinusoidal bases
              end                % with random phase
          end

          % Compute TV_LP coefficients (COVARIANCE type formulation, no window)

          sigL = length(sig);
          S = zeros(sigL,lpcRdr*(polRdr+1)); 	 % matrix for weighted past samples
	  s = sig(:);                            % Signal segment column
	  for k = 0:polRdr
	      basis(k+1) = 
	      for p = 1:lpcRdr                   % n-shift curFrm and ...
%		  temp = [prvFrm(end-p+1:end);curFrm(1:end-p)]; % insert past samples
		  temp = [zeros(p,1);sigL(1:end-p)]; % shifted signal with zero past
		  S(:,k*lpcRdr+p) = temp.*basis; % to get s[n-p].*bases_k[n]
	      end
	  end

	  tvLPwts = S\s;                         % Solve the Least Squares problem
	  tvLPres = s-S*tvLPwts;                 % least squares residual error

%	  EtvLP(cnt1) = tvLPres'*tvLPres;
% construct LP coefficient polynomial from the optimum weights

	  tvLPC = zeros(sigL,lpcRdr);
	  for p=1:lpcRdr
	      temp=zeros(sigL,1);
	      for k=0:polRdr
		  temp=temp + tvLPwts(k*lpcRdr+p)*basis(k+1);  % {sum( a_pk * n^k)}_k
	      end
	      tvLPC(:,p)=temp;
	  end
  return

%-----------------------------------------------------------------------