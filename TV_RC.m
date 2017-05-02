%------------------------------------------------------------
% Time-varying linear prediction (TVLP) analysis of speech
% using Lattice formulation and then extending it to 
% log-area ratios, to ensure stability of the filter
% parameters. The approach is based on the paper:
% "Autoregressive models with time-dependent log area-ratios"
% by Grenier, IEEE Tr ASSP, 1988.
% 
% Time-varying RC (K_i=1:p) parameters are further converted
% to  (i) TV_area function (+ve valued)
%     (ii) TV_predictor coefficients (A_k)
%     (iii) TV_LSF coefficients
% All these parametric solutions are better than taking the
% reverse route of starting from A_k, because the synthesis 
% filter stability is guaranteed. The TV solutions will not
% be same for different formulations of TV_LP or TV_Lattice 
% or TV_AreaF.
%
% Important experimental parameters are the (a) type of basis
% functions (polynomial basis is used) (b) number of basis 
% functions (polRdr) (c) order of prediction itself (lpcRdr).
% These parameters are dependent on the nature of the signal,
% sampling frequency and duration of the signal being modeled.
%
% The least squares formulation is similar to the covariance
% formulaton of linear prediction, hence no explicit window
% function is used. 
%
% All the TVLP parameters can then be resampled arbitrarily,
% because the parametric basis representation is available.
%
%  % the LP spectrum is computed for the spectrogram. The error
%  % signal power acts as a scale factor for the whole window,
%  % for comparison with the short-time signal spectrogram.
%----------------------------------------------------------------------------------------------
  function [RCpls RCmns RCavg resPls1 resMns1 resPls2 resMns2] = TV_RC(sig,lpcRdr,polRdr,type);
%----------------------------------------------------------------------------------------------

      sigL = length(sig);
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
              end              % with random phase
          end
      end
      
      for i = 0:polRdr,     % cross-product of bases for Phi matrix
          for j = 0:polRdr,
              CrsBases(i+1,j+1,:) = bases(i+1,:).*bases(j+1,:);
          end
      end
      
      resPls1 = sig;     % init Epsilon(t) for lattice forward prediction error
      resMns1 = sig;     % init Eeta(t) for lattice backward prediction error
      
      resPls2 = sig;     % similar inits for Burg-avg solution:
      resMns2 = sig;     % Kavg instead of Kpls and Kmns
      
      for p = 1:lpcRdr,     % successive lattice stages
          for i = 0:polRdr, % build 'Xi' and 'Hech' matrices
              Xim1(i+1,:) = bases(i+1,:).*resPls1'; % bases multiplied input
              Him1(i+1,:) = bases(i+1,:).*resMns1';
          end
          SQresPls1 = resPls1.^2;   % square of residue input for Phi matrices
          SQresMns1 = resMns1.^2;
          for i = 0:polRdr,
              for j = 0:polRdr,
                  for m = 1:sigL,
                      temp1(m) = CrsBases(i+1,j+1,m)*SQresPls1(m);
                      temp2(m) = CrsBases(i+1,j+1,m)*SQresMns1(m);
                  end
                  PhiPls(i+1,j+1) = sum(temp2(p+1:sigL));
                  PhiMns(i+1,j+1) = sum(temp1(p+1:sigL));
              end
              PsyPls(i+1) = resPls1(p+1:sigL)'*Him1(i+1,p:sigL-1)';
              PsyMns(i+1) = resMns1(p:sigL-1)'*Xim1(i+1,p+1:sigL)';
          end
          Kpls(p,:) = PhiPls\(-PsyPls'); % set up matrix equations & solve
          Kmns(p,:) = PhiMns\(-PsyMns'); % opt bases weights: Kpls_ij, Kmns_ij

          RCpls(p,:) = zeros(1,sigL);    % init for plus/minus solutions of K_ij
          RCmns(p,:) = zeros(1,sigL);
          for j = 0:polRdr,         % compute RC contours for Kpls, Kmns
              RCpls(p,:) = RCpls(p,:)+Kpls(p,j+1)*bases(j+1,:);
              RCmns(p,:) = RCmns(p,:)+Kmns(p,j+1)*bases(j+1,:);
          end

          ndx_pls = (abs(RCpls(p,:))>=1);  % chk for instability
          ndx_mns = (abs(RCmns(p,:))>=1);
          if sum(ndx_pls) > 0,            % of forward prediction
              ratio = 100*sum(ndx_pls)/sigL;
              fprintf('%5.0f %% unstable RCpls at p = %3.0f \n',ratio,p);
              for m = 1:sigL,
                  if ndx_pls(m) == 1,
                      RCpls(p,m) = sign(RCpls(p,m))*0.99; % correct instability
                  end
              end
          end
          if sum(ndx_mns) > 0,            % of backward prediction
              ratio = 100*sum(ndx_mns)/sigL;
              fprintf('%5.0f %% unstable RCmns at p = %3.0f \n',ratio,p);
              for m = 1:sigL,
                  if ndx_mns(m) == 1,
                      RCmns(p,m) = sign(RCmns(p,m))*0.99; % correct instability
                  end
              end
          end
                        % compute pth stage output for TV_lattice and
                        % then assign it as input to next stage iteration

          temp1 = resPls1 + [0 ; RCpls(p,1:sigL-1)'.* resMns1(1:sigL-1)];
          resMns1 = [0 ; resMns1(1:sigL-1)] + resPls1.*RCmns(p,:)';
          resPls1 = temp1;
          
% Also, solve for Burg-avg formulation; same first stage input but different
% successive stage outputs, because of a single average K_ij solution, as below

          for i = 0:polRdr, % build Xi and Hech matrices
              Xim1(i+1,:) = bases(i+1,:).*resPls2'; % bases modified inputs
              Him1(i+1,:) = bases(i+1,:).*resMns2';
          end
          SQresPls2 = resPls2.^2;
          SQresMns2 = resMns2.^2;
          for i = 0:polRdr,
              for j = 0:polRdr,
                  for m = 1:sigL,
                      temp1(m) = CrsBases(i+1,j+1,m).*SQresPls2(m);
                      temp2(m) = CrsBases(i+1,j+1,m).*SQresMns2(m);
                  end
                  PhiPls(i+1,j+1) = sum(temp2(p+1:sigL));
                  PhiMns(i+1,j+1) = sum(temp1(p+1:sigL));
              end
              PsyPls(i+1) = resPls2(p+1:sigL)'*Him1(i+1,p:sigL-1)';
              PsyMns(i+1) = resMns2(p:sigL-1)'*Xim1(i+1,p+1:sigL)';
          end
          
          Phi = PhiPls+PhiMns;  % avg of forward & backward residue matrices
          Psy = PsyPls+PsyMns;  % for Kavg_ij
          Kavg(p,:) = Phi\(-Psy');
          RCavg(p,:) = zeros(1,sigL);
          for j = 0:polRdr,     % compute RC contour for Kavg
             RCavg(p,:) = RCavg(p,:)+Kavg(p,j+1)*bases(j+1,:);
          end
          
          ndx_avg = (abs(RCavg(p,:))>=1);  % chk for instability
          if sum(ndx_avg) > 0,            % of Burg avg predictor
              ratio = 100*sum(ndx_avg)/sigL;
              fprintf('%5.0f %% unstable RCavg at p = %3.0f \n',ratio,p);
              for m = 1:sigL,
                  if ndx_avg(m) == 1,
                      RCavg(p,m) = sign(RCavg(p,m))*0.99; % correct instability
                  end
              end
          end
                        % compute pth stage output for Burg-avg lattice
                        % and assign it as input to next stage
                        
          temp2 = resPls2 + [0 ; RCavg(p,1:sigL-1)'.*resMns2(1:sigL-1)];
          resMns2 = [0 ; resMns2(1:sigL-1)] + resPls2.*RCavg(p,:)';
          resPls2 = temp2;

% Compute log Area-Ratio coefft contours, by a MSE fit to the stable regions
% of RCavg at each stage-p
                                      
      end    % all of lpcRdr no. of stages are computed using both formulations

    return
%---------------------------------------------------------------
