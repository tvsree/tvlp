%---------------------------------------------------------------
% Function to fit a chosen bases to the given contours of
% LP parameter estimation; the least squares formulation is 
% used. A fixed dimension vector contour is input and the
% same set of bases are fitted to each of the vector component 
% contours. 
% The function returns optimum basis weights for each of the 
% vector component contours and mse measure of goodness of fit.
%-----------------------------------------------------------------------
  function [optWts fitdVec avgMSEdb] = contourFit(parVec,basisRdr,type);
%-----------------------------------------------------------------------

      [vecDim contL] = size(parVec);     % each row is a contour

% construct basis functions

      if ~isempty(strfind(type,'leg')),
          if (polRdr>5), disp('Legendre polRdr > 5!'); return; end
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
          randPhase = [0 rand(1,basisRdr)*2*pi];      
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
      
% each dimension contour is fitted separately, to suit different ranges of values

      for i = 1:vecDim,
          for j = 0:basisRdr,   % signal-bases correlation vector Psi
              Psi(j+1) = parVec(i,:)*bases(j+1,:)';
          end
          optWts(i,:) = Phi \ Psi';
          
          fitdVec(i,:) = zeros(1,contL);
          for j = 0:basisRdr,
              fitdVec(i,:) = fitdVec(i,:)+optWts(i,j+1)*bases(j+1,:);
          end
          errCont = parVec(i,:)-fitdVec(i,:);
          mse(i) = errCont*errCont';
          vecPwr(i) = parVec(i,:)*parVec(i,:)';
          mseDB(i) = 10*log10(mse(i)/vecPwr(i));
      end
      avgMSEdb = sum(mseDB)/vecDim;
  return
%----------------------------------------------------------
