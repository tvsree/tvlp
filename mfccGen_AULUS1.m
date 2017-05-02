%----------------------------------------------------------
% To compute LPC spectrum using Auto/Covar formulation
% from a given wav file and fixed order predictor.
% From the LPC vector, solve for LSF, RC, Area parametric
% representations for comparison. 
% The input signal file (optionally re-sampled to Fs=8K;
% signal envelope is computed for recognizing speech
% end points, i.e., remove background silence from further
% parametric analysis and the endpoints are saved in a
% text file. (this can be used for testing stage of ASR)
% Optionally for training, the endpoints are input from 
% a label file, corresponding to the wav file.
%
% By changing the window size and overlap parameters,
% we can obtain either WB or NB quasi-stationary (QS) 
% analysis for LP. In either case Auto/Covar formulations
% can be used.
%
% For TV_LP analysis only covar formulation is
% justified. TV_LP estimated parameter contours are 
% re-sampled as required by the subsequent DNN classifiers.
%
% Optionally MFCC generation has been added;
% both MFCC parameters and log FilterBank energy vectors
% are saved into respective files, in the same required
% for KALDI.
%-----------------------------------------------------------
  function mfccGen;
  
    global Fs sigL wndSiz wndShft;
    Fs = 8000; nfft = 512; filCnt = 0;  % default values, updated from wavfile

        spkrDirID = fopen(['/prosjekt/dBase_OLLO/spkrLst'],'r');% spkr Directory list
        while ~feof(spkrDirID),
            curSpkr = fgetl(spkrDirID); disp(curSpkr);
            if isempty(curSpkr), break; end
            wavLstID = fopen(['/prosjekt/dBase_OLLO/' curSpkr '/wavLst'],'r');

            while ~feof(wavLstID),
                wavFil = fgetl(wavLstID);       % read data file name
                if isempty(wavFil), break; end; % exit on a blank line
                filCnt = filCnt+1;                
                [sigInp,FsInp] = audioread(wavFil);  % FsInp is given by the .wav
%      sig = resample(sigInp,Fs,FsInp);         % resample data by Fs/FsInp
                Fs = FsInp;                     % no sub-sampling in this case
                sig = sigInp;
                sigL = length(sig);
            
%      [bgnSmp,endSmp] = spchNdpt(sig,Fs);      % get bgn-end points of speech
%      bgntim = bgnSmp/Fs; endtim = endSmp/Fs;
%      ndpFil = wavFil ; ndpFil(end-2:end) = 'ndp';
%      ndpFilID = fopen(['./QSexpt2_16K/test/' ndpFil],'w');
%      fprintf(ndpFilID,'%10.0f',bgnSmp,endSmp,Fs,sigL);
%      fclose(ndpFilID);

% read the label file corresponding to current wav file

                [path,name,ext] = fileparts(wavFil);
                
                lblFil = ['/prosjekt/dBase_OLLO/OLLO2.0_LABELS_FORCED_ALIGNMENT/' ...
                          curSpkr '/' name '.label'];
                lblID = fopen(lblFil,'r');
                fgetl(lblID);
                for i = 1:3,
                    text = fgetl(lblID);
                    bgnEnd(i,:) = sscanf(text,'%i', 2);
                end
                bgnSmp = floor(bgnEnd(1,1)*10^(-7)*Fs);    % Phnm-1 begin
                endSmp = floor(bgnEnd(3,2)*10^(-7)*Fs);    % Phnm-3 end
                midSmp1 = floor(bgnEnd(1,2)*10^(-7)*Fs);   % Phnm-2 begin, Phnm-1 end
                midSmp2 = floor(bgnEnd(2,2)*10^(-7)*Fs);   % Phnm-3 begin, Phnm-2 end
                fclose(lblID);
                
% do framewise LP analysis
                
                wndSiz = floor(Fs*0.030);         % analysis parameters
                ovrlp = floor(Fs*0.015);
                wndShft = wndSiz-ovrlp;
                
%                  win = window(@hamming,wndSiz);    % for autocor analysis
%                  difsig = sig - 0.9*[0 sig(1:sigL-1)']'; % pre-emphasized signal

                nbCeps = 13 ; nbFbank = 24;       % MFCC analysis
                [ceps earMag] = MFCC(sig(bgnSmp:endSmp),wndSiz,wndShft,Fs,nbCeps,nbFbank);
                
%                  lpcRdr = floor(Fs/1000) + 2;
%                  lpcRdr = 12;                     % less than 2+Fs/1000
%                  frmCnt = 0;
%                  for n = bgnSmp:wndShft:endSmp,
%                      if (n+wndSiz-1 > sigL), break; end;
%                      frmCnt = frmCnt+1;
%                      sigw = difsig(n:n+wndSiz-1).*win;
%                      sigSpc = 20*log10(abs(fft(sigw,nfft)));
%            
% Autocorrelation analysis
%  
%                      [Ak,Ep] = lpc(sigw,lpcRdr);   % windowed signal for autocor-LP
%                      lpcSpc = -20*log10(abs(fft(Ak,nfft)))+10*log10(Ep*wndSiz);
%                      lsfAut = poly2lsf(Ak);        % convert LPCs to LSF [0:pi]
%                      lsfAutContr(:,frmCnt) = Fs*lsfAut/(2*pi);  % for plotting LSF in Hz Contours
%  
%                      rc = poly2rc(Ak);         % convert LPCs to reflection coeffs
%                      for i = 1:lpcRdr,         % chk for stability
%                          if abs(rc(i)) > 1,
%                              rc(i) = sign(rc(i))*0.99;
%                              fprintf('unstable Ak, at smp= %10.0f',n);
%                          end
%                      end
%                      rcContr(:,frmCnt) = rc;   % rc contour
%                      
%                      lar = rc2lar(rc);         % convert RCs to log-area-ratios (-ve also)
%                      larContr(:,frmCnt) = lar;
%                      
%                      areaRatio = exp(lar);     % convert to +ve valued area ratios
%                      area = cumprod(areaRatio);
%                      areaContr(:,frmCnt) = [1 ; area];

% Covariance analysis

%            [AkCov,EpCov] = arcov(difsig(n:n+wndSiz-1),lpcRdr);
%            covSpc = -20*log10(abs(fft(AkCov,nfft)))+10*log10(Ep*wndSiz);
%            lsfCov = lpc2lsf(AkCov,lpcRdr);
%            temp = size(lsfCov,2);
%            if temp < lpcRdr, 
%                lsfCov = [lsfCov zeros(1,lpcRdr-temp)]';
%                fprintf('unstable AkCov at smp= %5.0f',n);
%            end
%            lsfCovContr(:,frmCnt) = Fs*lsfCov/(2*pi);

%            spcGrm(:,frmCnt) = sigSpc(1:nfft/2);  % signal spectrogram
%            autLPgrm(:,frmCnt) = lpcSpc(1:nfft/2);% Aug_LP spectra for TF plot
%            covLPgrm(:,frmCnt) = covSpc(1:nfft/2);% Cov_LP spectra for TF plot

%                end     % end of current file

%      figure;
%      plot(sigSpc(1:nfft/2)); hold on; axis tight;
%      plot(lpcSpc(1:nfft/2),'r');
%      plot(covSpc(1:nfft/2),'g');
%      lsfScale = round(lsfAut*(nfft-1)/(2*pi));
%      spcMin = ones(lpcRdr,1)*min(sigSpc(1:nfft/2));
%      stem(lsfScale,spcMin);
%      lsf2Scale = round(lsfCov*(nfft-1)/(2*pi))
%      stem(lsf2Scale,spcMin,'r');

% save LSF,RC, LAR, MFCC and FBenergy parameters into respective files, in KALDI format

                intSyLbl = str2num(name(7:9));      % syllable label 1:150
                if intSyLbl <= 70,                    % compute phone labels for KALDI expts
                    tmp1 = mod(intSyLbl,14);
                    if tmp1 == 0, tmp1 = 14; end
                    centPhnLbl = tmp1-1;              % center Consonant Lbl of VCV 0:13 &
                    surrPhnLbl = floor((intSyLbl-1)/14)+14;   % surrounding vowel label: 14-18
                else
                    tmp1 = intSyLbl-71;
                    tmp2 = floor(tmp1/8)+1;
                    centPhnLbl = tmp2+13;             % center vowel Lbl of CVC 14:23
                    surrPhnLbl = rem(tmp1,8);         % surrounding consonant label: 0-7
                end
      
%                  lsfFil = wavFil ; lsfFil(end-2:end) = 'lsf';
%                  lsfFilID = fopen(['./QSexpt_Kaldi/' lsfFil],'w');
%                  rcFil = wavFil ; rcFil(end-2:end) = 'rc '; rcFil(end) = '';
%                  rcFilID = fopen(['./QSexpt_Kaldi/' rcFil],'w');
%                  larFil = wavFil ; larFil(end-2:end) = 'lar';
%                  larFilID = fopen(['./QSexpt_Kaldi/' larFil],'w');
                
                mfcFil = [name '.mfc'];
                mfcFilID = fopen(['/prosjekt/tvs/QSexpt_MFCC/' mfcFil],'w');
                
                fbeFil = [name '.fbe'];
                fbeFilID = fopen(['/prosjekt/tvs/QSexpt_MFCC/' fbeFil],'w');
                
                frmCnt = size(ceps,2);
                for i = 1:frmCnt,
                    PhnLbl = surrPhnLbl;
                    if (bgnSmp+(i-1)*wndShft >= midSmp1) & ...
                       (bgnSmp+(i-1)*wndShft <= midSmp2),
                        PhnLbl = centPhnLbl;
                    end

%                      fprintf(lsfFilID,lsfFil);
%                      fprintf(lsfFilID,'%8.0f',filCnt-1,i-1);
%                      fprintf(lsfFilID,'%8.0f',lsfAutContr(:,i));
%                      fprintf(lsfFilID,'%8.0f \n',PhnLbl);
%  
%                      fprintf(rcFilID,rcFil);
%                      fprintf(rcFilID,'%8.0f',filCnt-1,i-1);
%                      fprintf(rcFilID,'%8.3f',rcContr(:,i));
%                      fprintf(rcFilID,'%8.0f \n',PhnLbl);
%  
%                      fprintf(larFilID,larFil);
%                      fprintf(larFilID,'%8.0f',filCnt-1,i-1);
%                      fprintf(larFilID,'%8.3f',larContr(:,i));
%                      fprintf(larFilID,'%8.0f \n',PhnLbl);

                    fprintf(mfcFilID,mfcFil);
                    fprintf(mfcFilID,'%8.0f',filCnt-1,i-1);
                    fprintf(mfcFilID,'%11.3f',ceps(:,i));
                    fprintf(mfcFilID,'%8.0f \n',PhnLbl);

                    fprintf(fbeFilID,fbeFil);
                    fprintf(fbeFilID,'%8.0f',filCnt-1,i-1);
                    fprintf(fbeFilID,'%11.3f',earMag(:,i));
                    fprintf(fbeFilID,'%8.0f \n',PhnLbl);

                end
%                  fclose(larFilID);
%                  fclose(lsfFilID);
%                  fclose(rcFilID);
                  fclose(mfcFilID);
                  fclose(fbeFilID);

% Show TF plots of signal and LP parameters

%        time = [1:wndShft:sigL]/Fs;    % time,freq scales for TFR
%        freq = Fs*[0:nfft/2-1]/nfft;
%        timeLcl = [bgnSmp:wndShft:endSmp]/Fs; % silence removed region
%      
%        figure;
%        subplot(2,1,1); plot([1:sigL]/Fs,sig); axis tight; % axis in time units
%        subplot(2,1,2); stem([bgnSmp endSmp],[1 1],'fill');
%                        axis([1 sigL 0 1]);
%        figure;                               % plot signal spectrogram
%        surf(timeLcl,freq,spcGrm,'EdgeColor','none');
%        axis([time(1) time(end) freq(1) freq(end)]);
%        axis xy; colormap(jet); view(0,90);
%        xlabel('TIME (secs)'); ylabel('FREQUENCY (Hz)');
%        title('NarrowBand spectrogram');
%      
%        figure;                               % plot aut-LP spectrogram + LSFs
%        surf(timeLcl,freq,autLPgrm,'EdgeColor','none');
%        axis([time(1) time(end) freq(1) freq(end)]);
%        axis xy; colormap(jet); view(0,90);
%        xlabel('TIME (secs)'); ylabel('FREQUENCY (Hz)');
%        title('NB Autocor-LP spectrogram');
%        hold on;
%        for i = 1:lpcRdr,
%            plot(timeLcl,lsfAutContr(i,:),'m');
%        end
%        hold off;
%  
%        figure;                               % plot cov-LP spectrogram + LSFs
%        surf(timeLcl,freq,covLPgrm,'EdgeColor','none');
%        axis([time(1) time(end) freq(1) freq(end)]);
%        axis xy; colormap(jet); view(0,90);
%        xlabel('TIME (secs)'); ylabel('FREQUENCY (Hz)');
%        title('NB Covariance-LP spectrogram');
%        hold on;
%        for i = 1:lpcRdr,
%            plot(timeLcl,lsfCovContr(i,:),'m');
%        end
%        hold off;
     
            end         % end of 150 utterance files of a speaker
            fclose(wavLstID);
        end             % end of 50 speakers' list (5 dialects)
        fclose(spkrDirID);        
%    end                 % end of dialects
    
  return
%---------------------------------------------
