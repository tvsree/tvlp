%------------------------------------------------------------
% Generate time-varying parameters: TV_LPCs and TV_RCs
% using selected basis functions, over a longer time window.
% Two different formulations, TV_LP and TV_Lattice are used
% respectively. TV_LSFs are derived through LP->LSF or
% RC->LP->LSF transformation. The LP or RC bases weights
% are saved into respective files for KALDI based recognition
% experiments, to compare with quasi-stationary pars. Bases
% weights for TV_LSFs are also computed through post-fitting.
%
% OLLO dBase of CVC/VCV syllables and their label files;
% Fs = 16K is not modified; instead a lower lpcRdr = 12 is
% used to test.
% This code is run under the ../2016_NTNU/ directory
%------------------------------------------------------------
  function TVparGen;
%---------------------

    filCnt = 0; Fs = 8000; nfft = 512;

    dlctID = fopen('./OLLO/dialects');
    while ~feof(dlctID),
        curDlct = fgetl(dlctID);
        if isempty(curDlct), break; end
        spkrDirID = fopen(['./OLLO/' curDlct '/spkrLst'],'r');% spkr Directory list

        while ~feof(spkrDirID),
            curSpkr = fgetl(spkrDirID);
            if isempty(curSpkr), break; end
            wavLstID = fopen(['./OLLO/' curDlct '/' curSpkr '/wavLst'],'r');
        
            while ~feof(wavLstID),
                wavFil = fgetl(wavLstID);       % read data file name
                if isempty(wavFil), break; end; % exit on a blank line
                wavFilpath = ['./OLLO/' curDlct '/' curSpkr '/' wavFil];
                filCnt = filCnt+1;
                [sigInp,FsInp] = audioread(wavFilpath);  % FsInp is given by the .wav
                %      sig = resample(sigInp,Fs,FsInp);  % resample data by Fs/FsInp
                Fs = FsInp;                     % no sub-sampling in this case
                sig = sigInp;
                sigL = length(sig);

%        [bgnSmp,endSmp] = spchNdpt(sig,Fs); % get bgn-end points of speech
%        bgntim = bgnSmp/Fs; endtim = endSmp/Fs;
%        ndpFil = wavFil ; ndpFil(end-2:end) = 'ndp';
%        outFilID = fopen(['./QSexpt0_16K/' ndpFil],'w');
%        fprintf(outFilID,'%10.0f',bgnSmp,endSmp,Fs,sigL);
%        fclose(outFilID);

% read the label file corresponding to current wav file

                lblFil = [wavFil(1:end-3) 'label'];
                lblID = fopen(['./OLLO/OLLO2.0_LABELS_FORCED_ALIGNMENT/' ...
                               curSpkr '/' lblFil],'r');
                fgetl(lblID);
                for i = 1:3,
                    text = fgetl(lblID);
                    bgnEnd(i,:) = sscanf(text,'%i', 2);
                end
                bgnSmp = floor(bgnEnd(1,1)*10^(-7)*Fs);    % Phnm-1 begin
                endSmp = floor(bgnEnd(3,2)*10^(-7)*Fs);    % Phnm-3 end
                midSmp1 = floor(bgnEnd(1,2)*10^(-7)*Fs);   % Phnm-2 begin, Phnm-1 end
                midSmp2 = floor(bgnEnd(2,2)*10^(-7)*Fs);   % Phnm-3 begin, Phnm-2 end
      
                wndSiz = floor(Fs*0.030);         % analysis parameters
                ovrlp = floor(Fs*0.015);
                wndShft = wndSiz-ovrlp;

                %      lpcRdr = floor(Fs/1000) + 2;
                lpcRdr = 12;                      % less than (Fs/1000 + 2)
                wndSiz = floor(Fs*0.030);         % analysis parameters
                wndShft = floor(Fs*0.005);
                win = window(@hamming,wndSiz);    % for autocor analysis
                polRdr = 10;
      
                difsig = sig - 0.9*[0 sig(1:sigL-1)']'; % pre-emphasized difsig
     
                frmCnt1 = 0;
                for n = bgnSmp:wndShft:endSmp,
                    if (n+wndSiz-1 > sigL), break; end;
                    frmCnt1 = frmCnt1+1;
                    sigw = difsig(n:n+wndSiz-1).*win; % preemphasized windowed signal
                    sigSpc = 20*log10(abs(fft(sigw,nfft)));
                    spcGrm(:,frmCnt1) = sigSpc(1:nfft/2);  % signal spectrogram
          
% Autocorrelation analysis

                    [Ak,Ep] = lpc(sigw,lpcRdr);   % windowed signal for autocor-LP
                    lpcSpc = -20*log10(abs(fft(Ak,nfft)))+10*log10(Ep*wndSiz);
                    autLPgrm(:,frmCnt1) = lpcSpc(1:nfft/2);% Aut_LP spectra for TF plot
                    lsfAut = poly2lsf(Ak);        % convert LPCs to LSF [0:pi]
                    lsfAutContr(:,frmCnt1) = Fs*lsfAut/(2*pi);  % for plotting LSF Contours

                    rc = poly2rc(Ak);         % convert LPCs to reflection coeffs
                    for i = 1:lpcRdr,         % chk for stability
                        if abs(rc(i)) > 1,
                            fprintf('unstable Ak frame at %10.0f',n);
                            rc(i) = sign(rc(i))*0.99; % instability correction
                        end
                    end
%                    rcContr(:,frmCnt1) = rc;   % rc contour
%
%                    lar = rc2lar(rc);         % convert RCs to log area-ratios (-ve also)
%                    larContr(:,frmCnt1) = lar;
%                      areaRatio = exp(lar);  % convert to area ratios, +ve valued
%                      area = cumprod(areaRatio);
%                      areaContr(:,frmCnt1) = [1 ; area];
                end
                [fitdContr avgMSEdb] = contourFit(lsfAutContr,polRdr,'sin');

% TV_LPC analysis; whole syllable is analysed using a TVLP model

                [TVAk TVLPres] = TV_LPC(difsig(bgnSmp:endSmp),lpcRdr,polRdr,'sin');
                % choose among the bases: "pol"ynomial, "sin"e or "Leg"endre
                frmCnt2 = 0;
                for m = 1:wndShft:endSmp-bgnSmp+1, % sampling of continuous LPC contour
                    frmCnt2 = frmCnt2+1;
                    Ak_smp = TVAk(:,m);
                    resEnlocal = TVLPres(m:m+wndShft-1)'*TVLPres(m:m+wndShft-1);
                    TVlpcSpc = -20*log10(abs(fft(Ak_smp,nfft)))+10*log10(resEnlocal);
                    TVAkgrm(:,frmCnt2) = TVlpcSpc(1:nfft/2); % TVAk spectra for TF plot
                    
                    RCsmp = poly2rc(Ak_smp);       % convert LPCs to RCs and
                    for i = 1:lpcRdr,              % chk for stability
                        if abs(RCsmp(i)) > 1,
                            fprintf('unstable Ak_smp; m= %8.0f, i = %3.0f \n',m,i);
                            RCsmp(i) = sign(RCsmp(i))*0.99; % instability corrected
                        end
                    end
                    % rcContr(:,cnt2) = RCsmp;     % corrected RC contour
                    Ak_smp = rc2poly(RCsmp);       % recompute stable LPCs
                    lsf = lpc2lsf(Ak_smp,lpcRdr);  % convert LPCs to LSF [0:pi]
                    temp = size(lsf,2);
                    if temp < lpcRdr, 
                        lsf = [lsf zeros(1,lpcRdr-temp)]';
                        fprintf('unstable LSF at m= %5.0f \n',m);
                    end
                    TVLSF(:,frmCnt2) = Fs*lsf/(2*pi);  % convert LSF to Hz for plotting
                end                
                
% TV_Lattice analysis; whole syllable is analysed using the TV_Lattice model

                [RCpls RCmns RCavg resPls1 resMns1 resPls2 resMns2] = TV_Lattice(difsig(bgnSmp:endSmp),lpcRdr,polRdr,'sin');
                % choose among the bases: "pol"ynomial, "sin"e or "Leg"endre

                frmCnt3 = 0;
                for m = 1:wndShft:endSmp-bgnSmp+1,           % sampling of continuous contour
                    frmCnt3 = frmCnt3+1;
                    TVlpc(:,frmCnt3) = rc2poly(RCavg(:,m));  % convert TVRC at 'm' to LPC
                    resEnlocal = resPls1(m:m+wndShft-1)'*resPls1(m:m+wndShft-1);
                    TVlpcSpc = -20*log10(abs(fft(TVlpc(:,frmCnt3),nfft)))+10*log10(resEnlocal);
                    TVLPgrm(:,frmCnt3) = TVlpcSpc(1:nfft/2); % TVLP spectra for TF plot
                    
                    LSFpls = poly2lsf(TVlpc(:,frmCnt3));      % convert LPCs to LSF [0:pi]
                    TVLSFpls(:,frmCnt3) = Fs*LSFpls/(2*pi);  % for plotting LSF Contours
                    lar = rc2lar(RCavg(:,m));                % log area-ratio from RC contour
                    TVarea(:,frmCnt3) = cumprod(exp(lar));   % area function contour
                end

% Show TF plots of signal and LSF parameters, both QS and TV

                time = [1:wndShft:sigL]/Fs;    % time,freq scales for TFR
                freq = Fs*[0:nfft/2-1]/nfft;
                timeLcl = [bgnSmp:wndShft:endSmp]/Fs; % silence removed region
    
                figure;
                subplot(2,1,1); plot([1:sigL]/Fs,sig); axis tight; % axis in time units
                subplot(2,1,2); stem([bgnSmp endSmp],[1 1],'fill'); axis([1 sigL 0 1]);
%                  subplot(3,1,3); axis([time(1) time(end) freq(1) freq(end)]);
%                  hold on;
%                  for i = 1:lpcRdr,
%                      plot(timeLcl,TVLSFmns(i,:),'m'); 
%                  end
%                  hold off;
      
                figure;                               % plot signal spectrogram
                surf(timeLcl,freq,spcGrm,'EdgeColor','none');
                axis([time(1) time(end) freq(1) freq(end)]);
                axis xy; colormap(jet); view(0,90);
                xlabel('TIME (secs)'); ylabel('FREQUENCY (Hz)');
                title('NB-Spectrogram');
    
                figure;                               % plot QSaut-LP spectrogram
                surf(timeLcl,freq,autLPgrm,'EdgeColor','none');
                axis([time(1) time(end) freq(1) freq(end)]);
                axis xy; colormap(jet); view(0,90);
                xlabel('TIME (secs)'); ylabel('FREQUENCY (Hz)');
                title('QSaut-LP spectrogram + QS-LSF');
      
                hold on;                              % plot QS_LSF-contours
                for i = 1:lpcRdr,
                    plot(timeLcl,lsfAutContr(i,:),'m');
                    plot(timeLcl,fitdContr(i,:),'g')
                end
                hold off;
      
                figure;                               % plot TV_LPC spectrogram
                surf(timeLcl,freq,TVAkgrm,'EdgeColor','none');
                axis([time(1) time(end) freq(1) freq(end)]);
                axis xy; colormap(jet); view(0,90);
                xlabel('TIME (secs)'); ylabel('FREQUENCY (Hz)');
                title('TVLPC spectrogram + TVLSF');
      
                hold on;                              % plot TVLSF-contours
                for i = 1:lpcRdr,
                    plot(timeLcl,TVLSF(i,:),'m');
                end
                hold off;
                
                figure;                               % plot TV_RC spectrogram
                surf(timeLcl,freq,TVLPgrm,'EdgeColor','none');
                axis([time(1) time(end) freq(1) freq(end)]);
                axis xy; colormap(jet); view(0,90);
                xlabel('TIME (secs)'); ylabel('FREQUENCY (Hz)');
                title('TVRC spectrogram + TVLSFpls');
      
                hold on;                              % plot TVLSF-contours
                for i = 1:lpcRdr,
                    plot(timeLcl,TVLSFpls(i,:),'m');
                end
                hold off;

%                  figure;                               % plot QSautLP_Area spectrogram
%                  surf(timeLcl,[1:lpcRdr],larContr,'EdgeColor','none');
%                  axis([time(1) time(end) 0 lpcRdr+1]);
%                  axis xy; colormap(jet); view(0,90);
%                  xlabel('TIME (secs)'); ylabel('LAR: glottis->lips posn');
%                  title('QSaut-LP LAR spectrogram');
%        
%                  filCnt = 1; intLabel = 0;
%                  larFil = wavFil ; larFil(end-2:end) = 'lar'; % save LAR_Fn data
%                  outFilID = fopen(['./QSexpt0_16K/' larFil],'w');
%                  for i = 1:frmCnt,
%                      fprintf(outFilID,larFil);
%                      fprintf(outFilID,'%8.0f',filCnt-1,i-1);
%                      fprintf(outFilID,'%8.3f',larContr(:,i));
%                      fprintf(outFilID,'%8.0f \n',intLabel);
%                  end
%                  fclose(outFilID);
      
            end
        end
    end
  return
%---------------------------------------------------------------