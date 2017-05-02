clc;clear all;close all;

datafil='TVSkanDgts8K.wav';


[signal samplingRate]=audioread(datafil);signal=signal(:);
signal(2:end)=signal(2:end)-0.97*signal(1:end-1);

%%
windowSize=400*1e-3*samplingRate;
nFrames=floor(length(signal)/windowSize);

parameters.samplingRate=samplingRate;
parameters.nPoles=10;
parameters.ord=4;

prevFrame=zeros(windowSize,1);
for frameID=1:nFrames

    startPoint=(frameID-1)*windowSize+1;
    endPoint=startPoint+windowSize-1;
    signalFrame=signal(startPoint:endPoint);
    if(frameID>1)
        prevFrame=signal(startPoint-parameters.nPoles:1:startPoint-1);
    end
    
    %% For Least Squares TVLPC 
    [tvLPCParameters tvLPCResidual tvLPCCoefficients]=tvLPC(signalFrame,parameters,prevFrame,'LS');
    LS.tvLPCParameters=tvLPCParameters;
    LS.tvLPCResidual=tvLPCResidual;
    LS.tvLPCCoefficients=tvLPCCoefficients;
    
    %% For Sparse TVLPC
    [tvLPCParameters tvLPCResidual tvLPCCoefficients]=tvLPC(signalFrame,parameters,prevFrame,'SP');    
    SP.tvLPCParameters=tvLPCParameters;
    SP.tvLPCResidual=tvLPCResidual;
    SP.tvLPCCoefficients=tvLPCCoefficients;    
    

end