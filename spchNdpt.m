%-------------------------------------------------
% Compute short-time log-energy envelope 
% and the begin and end points of VCV/CVC speech
% segment, distinguishing from background silence
%-------------------------------------------------
  function [bgnSmp,endSmp] = spchNdpt(sig,Fs);
%----------------------------------------------
    sigL = length(sig);
    wndSiz = round(0.02*Fs); 	% 20 mSec window width
    wndShft = round(0.01*Fs);	% 10 mSec window shift
    minSeg = 0.1*Fs/wndShft;
    frmCnt = 0;
    for n = 1:wndShft:sigL-wndSiz,
	stPwr = (sig(n:n+wndSiz-1)'*sig(n:n+wndSiz-1))/wndSiz;
	frmCnt = frmCnt+1;
	Env(frmCnt) = stPwr;
    end
    maxPwr = max(Env);            % relative to max power ...
    EnvdB = 10*log10(Env/maxPwr); % usual speech power variation dB
    
    minPwrdB = min(EnvdB);        % min background power
    thr = max([minPwrdB+10 -40]); % for clean speech recording
    frmNdx = find(EnvdB > thr);   % locate speech seg > -40dB
    nSeg = 1;
    bgnSeg(nSeg) = frmNdx(1);
    for i = 2:size(frmNdx,2),
        if frmNdx(i)-frmNdx(i-1) > minSeg,
            width(nSeg) = frmNdx(i-1)-bgnSeg(nSeg)+1; % save segment widths
            nSeg = nSeg+1;
            bgnSeg(nSeg) = frmNdx(i);
        end
    end
    width(nSeg) = frmNdx(end)-bgnSeg(nSeg)+1;
    [mxWd ndx] = max(width);       % find longest segment
    bgnSmp = (bgnSeg(ndx)-1)*wndShft + 1; % frame posn -> samp index
    endSmp = bgnSmp-1 + width(ndx)*wndShft;
    
  return
%--------------------------------------------------------