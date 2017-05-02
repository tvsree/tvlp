%------------------------------------------------------------
% Time-varying linear prediction (TVLP) analysis of speech
% and display of TVLP spectrogram. TVLP solution formulation
% is based on polynomial basis functions of the predictor 
% coefficeints; the polynomial order and the predictor order
% is varied based on the signal window width and the signal
% sampling frequency, respectively.
%
% The least squares formulation is similar to the covariance
% formulaton of linear prediction, hence no explicit window
% function is used. 
%
% The optimum LP coefficients are sampled at every 'n' and
% the LP spectrum is computed for the LP_spectrogram. Error
% signal power at 'n' is scale factor for the whole window,
% for comparison with the short-time signal spectrogram.
%------------------------------------------------------------
  function tvLPspcGrm(filename);
%--------------------------------
      if (nargin < 1), filename = 'TVSkanDgts67_8K.wav'; end
      
      [sig, Fs, Nbits] = wavread(filename);
      winSiz = floor(Fs*0.2);
      skpSiz = floor(Fs*0.2);
      spcGrmSmp = Fs*0.001;
      sigL = length(sig);
      nfft = 512;
      lpcrdr = floor(Fs/1000) + 2;
      difsig = (sig - [0 sig(1:sigL-1)']');

      cnt1 = 0; cnt2 = 0;           % block count and frame count
      prvFrm = zeros(lpcrdr,1);
      resSig = [];
      for n = 1:skpSiz:sigL-winSiz,
	  cnt1 = cnt1+1;
%	  curFrm = difsig(n:n+winSiz-1);
%	  if n > 1, prvFrm = difsig(n-lpcrdr:n-1); end
  	  curFrm = sig(n:n+winSiz-1);
  	  if n > 1, prvFrm = sig(n-lpcrdr:n-1); end

% Compute TV_LP coefficients (COVARIANCE type formulation, no window)

	  nPoles = lpcrdr;                       % Number of poles
	  nOrd = 4;                              % Polynomial order
	  S = zeros(winSiz,nPoles*(nOrd+1)); 	 % matrix for past samples
	  s = curFrm(:);                         % Signal Frame

	  for k = 0:nOrd
	      basis = ([0:1:winSiz-1]/winSiz)'.^k;      % (n/N)^k for n=0:N-1
	      for p = 1:nPoles                   % n-shift curFrm and ...
		  temp = [prvFrm(end-p+1:end);curFrm(1:end-p)]; % insert past samples
		  S(:,k*nPoles+p) = temp.*basis; % to get s[n-p].*(n/N)^k
	      end
	  end

	  tvLPpars = S\s;                        % Solve the Least Squares problem
	  tvLPres = s-S*tvLPpars;                % LS residual error
	  EtvLP(cnt1) = tvLPres'*tvLPres;
	  resSig = [resSig tvLPres'];

% construct LP coefficient polynomial from the optimum weights

	  tvLPcoef=zeros(winSiz,nPoles);
	  for p=1:nPoles
	      temp=zeros(winSiz,1);
	      for k=0:nOrd
		  basis=([0:1:winSiz-1]/winSiz)'.^k;     % n^k
		  temp=temp+tvLPpars(k*nPoles+p)*basis;  % { sum( a_pk * n^k ) }_k
	      end
	      tvLPcoef(:,p)=temp;
	  end

% compute tvLP spectrogram for current block & append to previous block

	  for n = 1:spcGrmSmp:skpSiz,		 % ex: LP spectrum sampled @ 2.5ms = 20samp
	      cnt2 = cnt2+1;
	      tmp = tvLPres(n:n+spcGrmSmp-1);    % compute residue energy
	      EnLocal = tmp'*tmp;
	      Ak_smp = [1 -tvLPcoef(n,:)];       % LPcoef sampled at 'n'
	      tvLPspc = -20*log10(abs(fft(Ak_smp,nfft)))+10*log10(EnLocal);
	      tvLPspcGrm(:,cnt2) = tvLPspc(1:nfft/2);

% Also compute transformed parameters, RC, LSF, LAR from LP contour samples
	      
              RCsmp = poly2rc(Ak_smp);       % convert LPCs to RCs and
              for i = 1:lpcrdr,              % chk for stability
                  if abs(RCsmp(i)) > 1,
                      fprintf('unstable Ak_smp, at n= %8.0f \n',n);
                      RCsmp(i) = sign(RCsmp(i))*0.99; % instability corrected
                  end
              end
              rcContr(:,cnt2) = RCsmp;       % corrected RC contour
              Ak_smp = rc2poly(RCsmp);       % recompute stable LPCs
              
              lsfTVcov = poly2lsf(Ak_smp);   % compute LSF
              lsfTVcovContour(:,cnt2) = Fs*lsfTVcov/(2*pi);
              
              lar = rc2lar(RCsmp);           % compute log area-ratios (-ve also)
              larContr(:,cnt2) = lar;
              
              areaRatio = exp(lar);          % convert to area ratios, +ve valued
              area = cumprod(areaRatio);
              areaContr(:,cnt2) = [1 ; area];
	  end
      end
      
%       sigw = sig(n:n+winSiz-1);
%       sigSpc = 20*log10(abs(fft(sigw,nfft)));
%       [Akcov,Epcov] = arcov(sig(n:n+winSiz-1),lpcrdr);
%       covSpc = -20*log10(abs(fft(Akcov,nfft)))+10*log10(Epcov*winSiz);
% 
%       figure;
%       plot(sigSpc(1:nfft/2)); hold on; axis tight;
%       plot(covSpc(1:nfft/2),'r');
%       plot(tvLPspc(1:nfft/2),'g');
      
      figure;
      time = [1:spcGrmSmp:cnt1*skpSiz]/Fs;     % time, freq scales for TFR
      freq = Fs*[0:nfft/2-1]/nfft;             % display tvLP spectrogram
      surf(time,freq,tvLPspcGrm,'EdgeColor','none');
      axis xy; axis tight; colormap(jet); view(0,90);
      xlabel('TIME (secs)'); ylabel('FREQUENCY (Hz)');
      title('tvLPcov spectrogram');
      
      hold on;                                  % superpose LSF contours
      for i = 1:lpcrdr,
          plot(time,lsfTVcovContour(i,:),'m');
      end
      hold off;
      
      figure;                               % plot log area-ratio spectrogram
      time = [1:spcGrmSmp:cnt1*skpSiz]/Fs;  % time, Vtract posn scales for TFR
      freq = [1:lpcrdr];;
      surf(time,freq,larContr,'EdgeColor','none'); % display tvLAR spectrogram
      axis xy; axis tight; colormap(jet); view(0,90);
      xlabel('TIME (secs)'); ylabel('LAR: glottis-lips');
      title('tvLAR spectrogram');
      
      figure;                               % plot area-spectrogram
      time = [1:spcGrmSmp:cnt1*skpSiz]/Fs;  % time, Vtract posn scales for TFR
      freq = [0:lpcrdr];;
      surf(time,freq,areaContr,'EdgeColor','none'); % display tvLAR spectrogram
      axis xy; axis tight; colormap(jet); view(0,90);
      xlabel('TIME (secs)'); ylabel('Area: glottis-lips');
      title('tvArea spectrogram');
      
  return
      
%      figure;
%      subplot(3,1,1); plot(difsig(1:skpSiz*cnt1)); axis tight; title('source signal');
%      subplot(3,1,2); plot(resSig); title('tvLPres'); axis tight;

%      winSiz = floor(Fs*0.005);         % WBspec pars for covLP analysis
%      skpSiz = floor(Fs*0.0045);
%      covRes = zeros(1,sigL);           % covLP residue signal comparison
%      for n=lpcrdr+1:skpSiz:sigL-winSiz,
%          [Akcov,Epcov] = arcov(difsig(n:n+winSiz-1),lpcrdr);
%          for m = 0:skpSiz-1,
%              tmp = sig(n+m);
%              for k = 1:lpcrdr,
%                  tmp = tmp + difsig(n+m-k)*Akcov(k+1);
%              end
%              covRes(n+m) = tmp;
%          end
%      end
      
%      subplot(3,1,3); plot(covRes(1:length(resSig))); title('WBcovRes'); axis tight;

%-----------------------------------------------------------------------