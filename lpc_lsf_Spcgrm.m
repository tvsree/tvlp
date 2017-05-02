%-------------------------------------------------------
% To generate an LPC spectrum and LPC spectrogram
% from a given wav file. The LP analysis is performed
% using both autocorrelation and covariance analysis,
% with the same fixed order predictor, for comparison.
% The original spectrogram is also generated for the
% the same analysis parameters.
%
% By changing the window size and overlap parameters,
% we can obtain either a wide-band analysis or narrow-
% analysis for LP. In either case autocor or covar LP
% spectra are generated with proper residual signal
% scale factors to match with the original signal
% spectra.
%
% Spectrograms are plotted as three different figures
% for saving into eps files.
%-------------------------------------------------------
  function lpcSpcGrm(filename);
  if (nargin < 1), filename = 'TVSkanDgts67_8K.wav'; end
  
  [sig, Fs, Nbits] = wavread(filename);
  winSiz = floor(Fs*0.05);                % analysis parameters
  overlp = floor(Fs*0.045);
  skpSiz = winSiz-overlp;
  sigL = length(sig);
  nfft = 512;
  lpcrdr = floor(Fs/1000) + 2;
  win = window(@hamming,winSiz);          % for autocor analysis
  
  difsig = sig - 0.9*[0 sig(1:sigL-1)']'; % difference signal
  
  cnt = 0;
  for n = 1:skpSiz:sigL-winSiz,
      cnt = cnt+1;
      sigw = difsig(n:n+winSiz-1).*win;
      sigSpc = 20*log10(abs(fft(sigw,nfft)));

      [Ak,Ep] = lpc(sigw,lpcrdr);
      lpcSpc = -20*log10(abs(fft(Ak,nfft)))+10*log10(Ep*winSiz);

      lsfAut = poly2lsf(Ak);    % convert LPCs to LSF:[0:pi]
      lsfAutContour(:,cnt) = Fs*lsfAut/(2*pi);  % for plotting LSF Contour

      rc = poly2rc(Ak);         % convert LPCs to reflection coeffs, then
      lar = rc2lar(rc);         % convert RCs to log-area-ratios (-ve also)
      area = exp(lar);          % positive valued area function
      areaGrm(:,cnt) = lar;
      
      [Akcov,Epcov] = arcov(difsig(n:n+winSiz-1),lpcrdr);
      covSpc = -20*log10(abs(fft(Akcov,nfft)))+10*log10(Ep*winSiz);
      lsfCov = lpc2lsf(Akcov,lpcrdr);
      lsfCovContour(:,cnt) = Fs*lsfCov/(2*pi);

      autLPgrm(:,cnt) = lpcSpc(1:nfft/2); % save spectra for TF plot
      spcGrm(:,cnt) = sigSpc(1:nfft/2);
      covLPgrm(:,cnt) = covSpc(1:nfft/2);
      
%      figure;
%      plot(sigSpc(1:nfft/2)); hold on; axis tight;
%      plot(lpcSpc(1:nfft/2),'r');
%      plot(covSpc(1:nfft/2),'g');
%      lsfScale = round(lsfAut*(nfft-1)/(2*pi));
%      spcMin = ones(lpcrdr,1)*min(sigSpc(1:nfft/2));
%      stem(lsfScale,spcMin);
%      lsf2Scale = round(lsfCov*(nfft-1)/(2*pi))
%      stem(lsf2Scale,spcMin,'r');
      
  end

  time = [1:skpSiz:sigL-winSiz]/Fs;     % time,freq scales for TFR
  freq = Fs*[0:nfft/2-1]/nfft;
  
  figure;                        % plot signal spectrogram
  surf(time,freq,spcGrm,'EdgeColor','none');
  axis xy; axis tight; colormap(jet); view(0,90);
  xlabel('TIME (secs)'); ylabel('FREQUENCY (Hz)');
  title('NarrowBand spectrogram');
 
  figure;                        % plot aut-LP spectrogram + LSFs
  surf(time,freq,autLPgrm,'EdgeColor','none');
  axis xy; axis tight; colormap(jet); view(0,90);
  xlabel('TIME (secs)'); ylabel('FREQUENCY (Hz)');
  title('NB Autocor-LP spectrogram');
  hold on;
  for i = 1:lpcrdr,
      plot(time,lsfAutContour(i,:),'m');
  end
  hold off;

  figure;                        % plot cov-LP spectrogram + LSFs
  surf(time,freq,covLPgrm,'EdgeColor','none');
  axis xy; axis tight; colormap(jet); view(0,90);
  xlabel('TIME (secs)'); ylabel('FREQUENCY (Hz)');
  title('NB Covariance-LP spectrogram');
  hold on;
  for i = 1:lpcrdr,
      plot(time,lsfCovContour(i,:),'m');
  end
  hold off;

  figure;                        % plot aut-LP Log-area function
  surf(time,[1:lpcrdr],areaGrm,'EdgeColor','none');
  axis xy; axis tight; colormap(jet); view(0,90);
  xlabel('TIME (secs)'); ylabel('VocalTract-posn');
  title('NB autocor-LP log-area-ratio evolution');
  
%  spectrogram(sig,winSiz,overlp,nfft,Fs,'yaxis');

  return
%-----------------------------------------------------


%-----------------------------------------------------
% To compute the LP coefficients using the covariance
% formulation (non-stationary). Input signal is 
% unwindowed (implicit rectangular window) of length
% N samples. A p^th order optimum predictor is 
% computed. To determine the prediction error, a 
% Hamming windowed signal is passed through the LP
% filter and power is computed; this is done to compare
% the covarLP with the stationary LP formulation.
%-----------------------------------------------------
  function [Ak,Ep] = covarLP(x,p);
  
  N = length(x);
  for j = 1:p,
      for k = 1:p,
          tmp1 = 0;
	  for n = p+1:N,
	      tmp1 = tmp1+x(n-j)*x(n-k);
	  end
	  Rxx(j,k) = tmp1;	% cov matrix element
      end
      tmp2 = 0;
      for n = p+1:N,
          tmp2 = tmp2+x(n)*x(n-j);
      end
      Rx(j) = tmp2;	% right column vector
  end
  
  Ak = -Rxx\Rx';	% generalized matrix inversion

  w = hamming(N);	% Hamming window signal for
  xw = x.*w;		% residue signal computing
  res = zeros(N,1);	% to compare with autocor method
  for n = p+1:N,
      tmp = xw(n);
      for j = 1:p,
	  tmp = tmp+xw(n-j)*Ak(j);
      end
      res(n) = tmp;	% residue signal
  end
  Ep = res'*res;
  return
%-----------------------------------------------------
