  function test_vec;
      M = 10; N = 100; Rdr = 8;
      for m = 1:M,
          vec(m,:) = 2*m+sin(2*m*pi*[1:N]/N + pi*rand);
      end
      [optWts fitdVec avgMSEdb] = vecContFit(vec,Rdr,'sin');
      for m = 1:M,
          plot(vec(m,:),'b'); hold on;
          plot(fitdVec(m,:),'g');
      end
  end
%---------------------------------------------------------------
% Function to fit a chosen bases to the given contours of
% LP parameter estimation; the least squares formulation is 
% used. A fixed dimension vector contour is input and the
% same set of bases are fitted to each of the vector component 
% contours. 
% The function returns optimum basis weights for each of the 
% vector component contours and mse measure of goodness of fit.
%--------------------------------------------------------------------------
  function [optWts fitdVec avgMSEdb] = vecContFit(parVec,basisRdr,type);
%--------------------------------------------------------------------------

      [vecDim contL] = size(parVec);     % each row is a contour
      flag = 1;                 % for one time random phases

% construct basis functions

      if ~isempty(strfind(type,'leg')),
          if (basisRdr>5), disp('Legendre polRdr > 5!'); return; end
          legPol(1,:) = ones(1,contL);
          legPol(2,:) = [-contL/2:(contL/2-1)]*2/contL;
          legPol(3,:) = (3*legPol(2,:).^2 - 1)/2;
          legPol(4,:) = (5*legPol(2,:).^3 - 3*legPol(2,:))/2;
          legPol(5,:) = (35*legPol(2,:).^4 - 30*legPol(2,:).^2 + 3)/8;
          legPol(6,:) = (63*legPol(2,:).^5 - 70*legPol(2,:).^3 + 15*legPol(2,:))/8;
          
          for i = 1:basisRdr+1,
              bases(i,:) = legPol(i,:);
          end
      else
          if (flag),            % one time randomized phases
              randPhase = [0 rand(1,basisRdr)*2*pi];
              flag = 0;
          end
          for i = 0:basisRdr,   % polynomial bases upto basisRdr
              if ~isempty(strfind(type,'pol')), 
                  bases(i+1,:) = ([0:contL-1]/contL).^i; % polynomial bases
              end
              if ~isempty(strfind(type,'sin')), 
                  bases(i+1,:) = cos((i*2*pi*[0:contL-1]/contL) + randPhase(i+1)); % sinusoidal bases
              end               % with random phase
          end
      end
      for i = 0:basisRdr,       % cross-product of bases for Phi matrix
          for j = 0:basisRdr,
              Phi(i+1,j+1) = bases(i+1,:)*bases(j+1,:)';
          end
      end
      Phi_block = kron(eye(vecDim),Phi);
      
% all dimension contours are fitted jointly but separate weights for each 

      for i = 1:vecDim,
          tmp = (i-1)*(basisRdr+1);
          for j = 0:basisRdr,   % signal-bases correlation vector Psi
              Psi_tall(tmp+j+1) = parVec(i,:)*bases(j+1,:)';
          end
      end

      optWts_tall = Phi_block \ Psi_tall';
          
      for i = 1:vecDim,
          fitdVec(i,:) = zeros(1,contL);
          tmp = (i-1)*(basisRdr+1);
          for j = 0:basisRdr,
              fitdVec(i,:) = fitdVec(i,:)+optWts_tall(tmp+j+1)*bases(j+1,:);
          end
          errCont = parVec(i,:)-fitdVec(i,:);
          mse(i) = errCont*errCont';
          vecPwr(i) = parVec(i,:)*parVec(i,:)';
          mseDB(i) = 10*log10(mse(i)/vecPwr(i));
          optWts(i,:) = optWts_tall(tmp+1:tmp+basisRdr+1);
      end
      avgMSEdb = sum(mseDB)/vecDim;

      return
  end
%----------------------------------------------------------
