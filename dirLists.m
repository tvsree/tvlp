%----------------------------------------------
% Directory traversing and file management
% 1. Create file list of wav files by traversing 
%    speaker directories
% 2. Run a set of linux shell-scripts as 
%    required
%----------------------------------------------
    function dirLists;
%----------------------
        spkrLstID = fopen(['/prosjekt/dBase_OLLO/spkrLst'],'r');% spkr Directory list
        while ~feof(spkrLstID),
            curSpkr = fgetl(spkrLstID);
            if isempty(curSpkr), break; end
            curSpkrDir = ['/prosjekt/dBase_OLLO/' curSpkr '/'];
            eval(['!ls ' curSpkrDir '*_L*.wav > ' curSpkrDir 'wavLst']);
        end
    return
%-----------------------------------------------