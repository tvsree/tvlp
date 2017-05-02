function [filterWeights,fftFreqs] = fbankT(nbFilters,fftSize,samplingRate,minFreq,maxFreq)
% [fbTrMx,Freqs] =
%       fbankT(nbFilters,<fftSize,samplingRate,minFrequency,maxFrequency>)
% 
% Filter bank Transformation matrix  ( nbFilters x fftSize/2 )
% 
% defaults for optional parameters < >:
%   
%   fftSize      = 512
%   samplingRate = 8000
%   minFrequency = 0
%   maxFrequency = samplingRate/2
% 
%                                        -mijail. 20/march/2001
  
% FBANK Parameters: ********************
  
 if (nargin < 2) fftSize = 512; end;
 if (nargin < 3) samplingRate = 8000; end;
 if (nargin < 4) minFreq = 0; end;
 if (nargin < 5) maxFreq = samplingRate/2; end;


% PROCESSING : ********************************

 % Figure out the band edges.
 % Interesting frequencies are lineary spaced in Mel scale. 
 freqs=imel(linspace(mel(minFreq),mel(maxFreq),nbFilters+2));

 % Lower, center, and upper band edges are consecutive interesting freqs. 
 lower = freqs(1:nbFilters);
 center = freqs(2:nbFilters+1);
 upper = freqs(3:nbFilters+2);

 % Reserving memory for the transformation matrix
 filterWeights = zeros(nbFilters,fftSize/2);

 % Assuming a triangular weighting function.
 triangleHeight =ones(1,nbFilters);    % height is constant = 1
 % triangleHeight = 2./(upper-lower);  % weight is constant = 1

 % frequency bins
 fftFreqs = (0:fftSize/2-1)/fftSize*samplingRate;

 % Figure out each frequencies contribution
 for chan=1:nbFilters
	filterWeights(chan,:) =... 
  (fftFreqs > lower(chan) & fftFreqs <= center(chan)).* ...
   triangleHeight(chan).*(fftFreqs-lower(chan))/(center(chan)-lower(chan)) + ...
  (fftFreqs > center(chan) & fftFreqs < upper(chan)).* ...
   triangleHeight(chan).*(upper(chan)-fftFreqs)/(upper(chan)-center(chan));
end

% plot(fftFreqs,filterWeights');
% axis([lower(1) upper(nbFilters) 0 max(max(filterWeights))])
