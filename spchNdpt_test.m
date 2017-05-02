%-------------------
  function TestNdpt;
%-------------------
    Fs = 8000;
    path = './OLLO2.0_NO/S01F_NO/';
    wavFil = [path 'S01F_L140_V2_M1_N1_CS0.wav'];
    [sigInp,FsInp,Nbits] = wavread(wavFil);  % Fs is given by the .wav
%    sig = resample(sigInp,Fs,FsInp);    % resample data by Fs/FsInp
    Fs = FsInp; 
    sig = sigInp;
    [bgnSmp,endSmp] = spchNdpt(sig,Fs);   % get start & end points of speech
    sigL = length(sig);

    figure;
    subplot(3,1,1); plot([1:sigL]/Fs,sig); axis tight; % axis in time units
%    subplot(3,1,2); plot([1:frmCnt],EnvdB); axis tight;
    subplot(3,1,3); stem([bgnSmp endSmp],[1 1],'fill');
                      axis([1 sigL 0 1]);
  return
%--------------------------------------------------------------------