function [ceps earMag] = MFCC(Signal,windowSize,frameShift,samplingRate,nbCeps,nbFbank,fftSize,minFreq,maxFreq)
% MFCC - Mel frequency cepstrum coefficient analysis.
%
% [ceps,earMag] =
%      MFCC(Signal,<windowSize,frameShift,samplingRate,
%                   nbCeps,nbFbank,fftSize,minFreq,maxFreq>)
%
% defaults for optional parameters
%
%   windowSize   = 256
%   frameShift   = windowSize/2
%   samplingRate = 8000
%   nbCeps       = 12
%   nbFbank      = 40
%   fftSize      = 512
%   minFrequency = 0
%   maxFrequency = samplingRate/2
%
%                                       -mijail, 20/march/2001

% MFCC parameters: *******************

 if (nargin < 2) windowSize = 256; end;
 if (nargin < 3) frameShift = windowSize/2; end;
 if (nargin < 4) samplingRate = 8000; end;
 if (nargin < 5) nbCeps = 12; end;
 if (nargin < 6) nbFbank = 40; end;
 if (nargin < 7) fftSize = 512; end;
 if (nargin < 8) minFreq = 0; end;
 if (nargin < 9) maxFreq = samplingRate/2; end;

 a=-0.97;   % pre-emphasis filter coefficient


% PROCESSING  *********************

 % Filter the input Signal with the preemphasis filter.
 preEmphasized = filter([1 a], 1, Signal);

 % Windowing the signal
 windowedSignal= windowize(preEmphasized,windowSize,frameShift);

 % Weight the signal with a hamming window
 hamWindow = 0.54 - 0.46*cos(2*pi*(0:windowSize-1)'/windowSize);

 weightedSignal=hamWindow(:,ones(1,size(windowedSignal,2))).*windowedSignal;

 % Compute the Fourier Transform
 fftData=fft(weightedSignal,fftSize);

 % Take the amplitude of interesting coefficients
 fftMag=abs(fftData(1:fftSize/2,:));

 
 % Compute the log-FilterBank amplitudes
 earMag=log10(fbankT(nbFbank,fftSize,samplingRate,minFreq,maxFreq)*fftMag);

 % apply the Discrete Cosine Transform
 ceps= dct(earMag); ceps=ceps(1:nbCeps,:);


% % Alternately, we can compute the DCT transform as in the HTK toolkit.
% %
% % Figure out Discrete Cosine Transform.
% % dct(i,j) is a matrix  nbCeps x nbFbank.
% % The i,j component is given by
% %        dct(i,j)=sqrt(2/nbFbank)  *  cos( pi/nbFbank * i*(j+0.5) )
% % where we have assumed that i and j start at 0.
%
% dctMatrix = realsqrt(2./nbFbanks).*...
%    cos((0:nbCeps-1)'*((0:nbFbank-1)+0.5).*pi./nbFbank);
% ceps = dctMatrix * earMag
