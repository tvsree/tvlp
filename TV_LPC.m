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
  function [tvLPC tvLPres] = TV_LPC(sig,lpcRdr,polRdr,type);
%----------------------------------------------------------

      sigL = length(sig);

% construct basis functions

      if ~isempty(strfind(type,'leg')),
          if (polRdr>5), disp('Legendre polRdr > 5!'); return; end
          legPol(1,:) = ones(1,sigL);
          legPol(2,:) = [-sigL/2:(sigL/2-1)]*2/sigL;
          legPol(3,:) = (3*legPol(2,:).^2 - 1)/2;
          legPol(4,:) = (5*legPol(2,:).^3 - 3*legPol(2,:))/2;
          legPol(5,:) = (35*legPol(2,:).^4 - 30*legPol(2,:).^2 + 3)/8;
          legPol(6,:) = (63*legPol(2,:).^5 - 70*legPol(2,:).^3 + 15*legPol(2,:))/8;
          
          for i = 1:polRdr+1,
              bases(i,:) = legPol(i,:);
          end
      else
          randPhase = [0 rand(1,polRdr)*2*pi];      
          for i = 0:polRdr,     % polynomial bases upto polRdr
              if ~isempty(strfind(type,'pol')), 
                  bases(i+1,:) = ([0:sigL-1]/sigL).^i; % polynomial bases
              end
              if ~isempty(strfind(type,'sin')), 
                  bases(i+1,:) = cos((i*2*pi*[0:sigL-1]/sigL) + randPhase(i+1)); % sinusoidal bases
              end                % with random phase
          end
      end
            
% Compute TV_LP coefficients (COVARIANCE type formulation, no window)

      S = zeros(sigL,lpcRdr*(polRdr+1)); 	 % matrix for weighted past samples
      s = sig(:);                            % Signal segment column

      for k = 0:polRdr,
	  for p = 1:lpcRdr                   % n-shift curFrm and ...
	      temp = [zeros(p,1);sig(1:end-p)]; % shifted signal with zero past
	      S(:,k*lpcRdr+p) = temp.*bases(k+1,:)'; % to get s[n-p].*bases_k[n]
	  end
      end
      tvLPwts = S\s;                         % Solve the Least Squares problem
      tvLPres = s-S*tvLPwts;                 % least squares residual error

%	  EtvLP(cnt1) = tvLPres'*tvLPres;
% construct LP coefficient polynomial from the optimum weights

      tvLPC = zeros(lpcRdr+1,sigL);
      tvLPC(1,:) = ones(1,sigL);
      for p=1:lpcRdr
          temp=zeros(1,sigL);
          for k=0:polRdr
              temp=temp + tvLPwts(k*lpcRdr+p)*bases(k+1,:);  % {sum( a_pk * n^k)}_k
          end
	  tvLPC(p+1,:) = -temp;
      end
      
%            ndx_pls = (abs(RCpls(p,:))>=1);  % chk for instability
%            if sum(ndx_pls) > 0,            % of forward prediction
%                ratio = 100*sum(ndx_pls)/sigL;
%                fprintf('%5.0f %% unstable RCpls at p = %3.0f \n',ratio,p);
%                for m = 1:sigL,
%                    if ndx_pls(m) == 1,
%                        RCpls(p,m) = sign(RCpls(p,m))*0.99; % correct instability
%                    end
%                end
%            end

  return

%-----------------------------------------------------------------------