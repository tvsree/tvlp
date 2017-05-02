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
  function TVRCwtsGen;
%----------------------

    filCnt = 0; Fs = 8000; nfft = 512;

% select speech file from dBase

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
                [sigInp,FsInp] = audioread(wavFilpath); % FsInp is given by the .wav
                %      sig = resample(sigInp,Fs,FsInp); % resample data by Fs/FsInp
                Fs = FsInp;                             % no sub-sampling in this case
                sig = sigInp;
                sigL = length(sig);
      
                difsig = sig - 0.9*[0 sig(1:sigL-1)']'; % pre-emphasized difsig

%        [bgnSmp,endSmp] = spchNdpt(sig,Fs);            % get bgn-end points of speech
%        bgntim = bgnSmp/Fs; endtim = endSmp/Fs;
%        ndpFil = wavFil ; ndpFil(end-2:end) = 'ndp';
%        outFilID = fopen(['./QSexpt0_16K/' ndpFil],'w');
%        fprintf(outFilID,'%10.0f',bgnSmp,endSmp,Fs,sigL);
%        fclose(outFilID);

% read the label file corresponding to current wav file for KALDI expts

                [path,name,ext] = fileparts(wavFil);

                lblFil = [name '.label'];
                lblID = fopen(['./OLLO/OLLO2.0_LABELS_FORCED_ALIGNMENT/' ...
                               curSpkr '/' lblFil],'r');
%                lblID = ['/prosjekt/dBase_OLLO/OLLO2.0_LABELS_FORCED_ALIGNMENT/' ...
%                          curSpkr '/' lblFil];
                fgetl(lblID);
                for i = 1:3,
                    text = fgetl(lblID);
                    bgnEnd(i,:) = sscanf(text,'%i', 2);
                end
                bgnSmp = floor(bgnEnd(1,1)*10^(-7)*Fs);    % Phnm-1 begin
                endSmp = floor(bgnEnd(3,2)*10^(-7)*Fs);    % Phnm-3 end
                midSmp1 = floor(bgnEnd(1,2)*10^(-7)*Fs);   % Phnm-2 begin, Phnm-1 end
                midSmp2 = floor(bgnEnd(2,2)*10^(-7)*Fs);   % Phnm-3 begin, Phnm-2 end
                
                intSyLbl = str2num(name(7:9));        % syllable label 1:150
                if intSyLbl <= 70,                    % compute phone labels for KALDI expts
                    tmp1 = mod(intSyLbl,14);
                    if tmp1 == 0, tmp1 = 14; end
                    centPhnLbl = tmp1-1;              % center Consonant Lbl of VCV 0:13 &
                    surrPhnLbl = floor((intSyLbl-1)/14)+14; % surrounding vowel label: 14-18
                else
                    tmp1 = intSyLbl-71;
                    tmp2 = floor(tmp1/8)+1;
                    centPhnLbl = tmp2+13;             % center vowel Lbl of CVC 14:23
                    surrPhnLbl = rem(tmp1,8);         % surrounding consonant label: 0-7
                end
                
% start parameter estimation

                %      lpcRdr = floor(Fs/1000) + 2;
                lpcRdr = 12;                      % less than (Fs/1000 + 2)
                wndSiz = floor(Fs*0.030);         % analysis parameters
                wndShft = floor(Fs*0.015);
                win = window(@hamming,wndSiz);    % for autocor analysis

                polRdr = 10;
                bigWnd = wndSiz*11/2;           % signal size for TV analysis
                bigShft = wndSiz*5/2;           % independent of overlap in QS analysis
     
                frmCnt0 = 0;
                for n = bgnSmp:wndShft:endSmp-wndShft,
                    if (n+wndSiz-1 > sigL), break; end;
                    frmCnt0 = frmCnt0+1;
                    sigw = difsig(n:n+wndSiz-1).*win; % preemphasized windowed signal
                    sigSpc = 20*log10(abs(fft(sigw,nfft)));
                    spcGrm(:,frmCnt0) = sigSpc(1:nfft/2);  % signal spectrogram
          
% Autocorrelation analysis

                    [Ak,Ep] = lpc(sigw,lpcRdr);   % windowed signal for autocor-LP
                    lpcSpc = -20*log10(abs(fft(Ak,nfft)))+10*log10(Ep*wndSiz);
                    autLPgrm(:,frmCnt0) = lpcSpc(1:nfft/2);% Aut_LP spectra for TF plot
                    lsfAut = poly2lsf(Ak);        % convert LPCs to LSF [0:pi]
                    lsfAutContr(:,frmCnt0) = Fs*lsfAut/(2*pi);  % for plotting LSF Contours

                    rc = poly2rc(Ak);         % convert LPCs to reflection coeffs
                    for i = 1:lpcRdr,         % chk for stability
                        if abs(rc(i)) > 1,
                            fprintf('unstable Ak frame at %10.0f',n);
                            rc(i) = sign(rc(i))*0.99; % instability correction
                        end
                    end
                    rcContr(:,frmCnt0) = rc;   % rc contour

%                    lar = rc2lar(rc);         % convert RCs to log area-ratios (-ve also)
%                    larContr(:,frmCnt0) = lar;
%                      areaRatio = exp(lar);  % convert to area ratios, +ve valued
%                      area = cumprod(areaRatio);
%                      areaContr(:,frmCnt0) = [1 ; area];
                end
                
% TV_Lattice analysis; whole syllable is analysed using the TV_Lattice model
% open LPwts, LSFwts & other TVparameter files, to suit KALDI expts
%               outPath = '/prosjekt/tvs/TVRCexpt/'; for SIRKUS7

                outPath1 = './TVRCexpt/';
                tvlsfFil = [name '.tvlsf'];
                tvlsfFilID = fopen([outPath1 tvlsfFil],'w');
                tvrcFil = [name '.tvrc'];
                tvrcFilID = fopen([outPath1 tvrcFil],'w');
                %                  tvlarFilID = fopen([outPath name '.tvlar'],'w');
                outPath2 = './TVRCWts/';
                tvLwtsFil = [name '.tvLSFwts'];
                tvLwtsFilID = fopen([outPath2 tvLwtsFil],'w');
                tvRwtsFil = [name '.tvRCwts'];
                tvRwtsFilID = fopen([outPath2 tvRwtsFil],'w');
                %                  tvlarFilID = fopen([outPath2 name '.tvlar'],'w');

                frmCnt3 = 0; frmCnt5 = 0; prvCnt = 0;
                for n = bgnSmp:bigShft:endSmp,   % Big window TV_Lattice analysis
                    frmCnt3 = frmCnt3+1;
                    [Kpls Kmns Kavg RCpls RCmns RCavg resPls1 resMns1 resPls2 resMns2] = ...
                    TV_RCwts(difsig(n:n+bigWnd-1),lpcRdr,polRdr,'sin');
                    % choose among bases type: "pol"ynomial, "sin"e or "Leg"endre
                    
                    for m = 1:wndShft:bigShft,
                        frmCnt5 = frmCnt5+1;
                        TVAK(:,frmCnt5) = rc2poly(RCavg(:,m));   % convert TVRC at 'm' to LPC
                        resEnlocal = resPls2(m:m+wndShft-1)'*resPls2(m:m+wndShft-1);
                        TVAKSpc = -20*log10(abs(fft(TVAK(:,frmCnt5),nfft)))+10*log10(resEnlocal);
                        TVAKgrm(:,frmCnt5) = TVAKSpc(1:nfft/2);  % TVLP spectra for TF plot
                    
                        LSFpls = poly2lsf(TVAK(:,frmCnt5));      % convert LPCs to LSF [0:pi]
                        TVLSFpls(:,frmCnt5) = Fs*LSFpls/(2*pi);         % normalized LSF 0-1
                        
                        lar = rc2lar(RCavg(:,m));                % log area-ratio from RC contour
                        TVarea(:,frmCnt5) = cumprod(exp(lar));   % area function contour
                        
                        PhnLbl = surrPhnLbl;
                        if (n+m-1 >= midSmp1) & (n+m-1 <= midSmp2),  % get current window phone label
                           PhnLbl = centPhnLbl;
                        end
                        fprintf(tvlsfFilID,tvlsfFil);            % save TVLSF parameters
                        fprintf(tvlsfFilID,'%9.0f',filCnt-1,frmCnt5-1);
                        fprintf(tvlsfFilID,'%8.4f',LSFpls/pi);
                        fprintf(tvlsfFilID,'%8.0f \n',PhnLbl);

                        fprintf(tvrcFilID,tvrcFil);              % save TVRC parameters
                        fprintf(tvrcFilID,'%9.0f',filCnt-1,frmCnt5-1);
                        fprintf(tvrcFilID,'%8.3f',RCavg(:,m));
                        fprintf(tvrcFilID,'%8.0f \n',PhnLbl);

%                          fprintf(tvlarFilID,tvlarFil);
%                          fprintf(tvlarFilID,'%9.0f',filCnt-1,frmCnt5-1);
%                          fprintf(tvlarFilID,'%8.3f',larContr(:,frmCnt5));
%                          fprintf(tvlarFilID,'%8.0f \n',PhnLbl);

                    end
                
                    [LSFoptWts fitdContr avgMSEdb] = contourFit(TVLSFpls(:,prvCnt+1:frmCnt5),polRdr,'sin');
                    % choose among bases type: "pol"ynomial, "sin"e or "Leg"endre
                    prvCnt = frmCnt5;
                    
                    PhnLbl = surrPhnLbl;
                    if (n+bigWnd >= midSmp1) & (n <= midSmp2), ... % determine current phone label
                       PhnLbl = centPhnLbl;
                    end

                    fprintf(tvLwtsFilID,tvLwtsFil);
                    fprintf(tvLwtsFilID,'%9.0f',filCnt-1,frmCnt3-1);
                    fprintf(tvLwtsFilID,'%8.4f',LSFoptWts);
                    fprintf(tvLwtsFilID,'%8.0f \n',PhnLbl);
                    
                    fprintf(tvRwtsFilID,tvRwtsFil);
                    fprintf(tvRwtsFilID,'%9.0f',filCnt-1,frmCnt3-1);
                    fprintf(tvRwtsFilID,'%8.4f',Kavg);
                    fprintf(tvRwtsFilID,'%8.0f \n',PhnLbl);
                end
%                fclose(tvlarFilID);
                fclose(tvlsfFilID);
                fclose(tvrcFilID);
                fclose(tvLwtsFilID);
                fclose(tvRwtsFilID);
                
% Show TF plots of signal and LSF parameters, both QS and TV

                time = [1:wndShft:sigL]/Fs;    % time,freq scales for TFR
                freq = Fs*[0:nfft/2-1]/nfft;
                timeLcl = [bgnSmp:wndShft:endSmp-wndShft]/Fs; % silence removed region
    
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
      
                hold on;                              % superpose QS_LSF-contours
                for i = 1:lpcRdr,
                    plot(timeLcl,lsfAutContr(i,:),'m');
%                    plot(timeLcl,fitdContr(i,:),'g')
                end
                hold off;
      
%                  figure;                               % plot TV_LPC spectrogram
%                  surf(timeLcl,freq,TVAkgrm,'EdgeColor','none');
%                  axis([time(1) time(end) freq(1) freq(end)]);
%                  axis xy; colormap(jet); view(0,90);
%                  xlabel('TIME (secs)'); ylabel('FREQUENCY (Hz)');
%                  title('TVLPC spectrogram + TVLSF');
%        
%                  hold on;                              % plot TVLSF-contours
%                  for i = 1:lpcRdr,
%                      plot(timeLcl,TVLSF(i,:),'m');
%                  end
%                  hold off;
                
                figure;                               % plot TVRC->LP spectrogram
                surf(timeLcl,freq,TVAKgrm(:,1:frmCnt0),'EdgeColor','none');
                axis([time(1) time(end) freq(1) freq(end)]);
                axis xy; colormap(jet); view(0,90);
                xlabel('TIME (secs)'); ylabel('FREQUENCY (Hz)');
                title('TVRC->LP spectrogram -> TVLSFpls');
                
                hold on;                              % plot TVLSF-contours
                for i = 1:lpcRdr,
                    plot(timeLcl,TVLSFpls(i,1:frmCnt0),'m');
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
