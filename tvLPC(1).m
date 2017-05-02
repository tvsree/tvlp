%% Time Vaying Linear Prediction
% This function is used to estimate TVLP coefficients given the signal
% s[n]= { sum { sum( a_pk n^k s[n-p] ) }_p }_k + e[n]
% { sum () }_k denotes sum over k
% k: polynomial order p: number of poles
% Inputs: 
% currentFrame: s=signal(0:1:N-1)
% prevFrame: signal(-N:1:-1)
% parameters: It is structure containing the below variables
%   sampling rate: parameters.samplingRate (default:8000)
%   number of Poles: parameters.nPoles (default:samplingRate(kHz) +2)
%   polynomial order: parameters.ord (default: 3)
%
% Outputs:
% tvLPCParameters:  nPoles*(ord+1) vector
%    [a_10 a_20 ... a_p0 a_11 a_21 a_31 ... a_p1 ... ... a_1k a_2k ... a_pk]
% tvLPCResidual: vector of same size as speech segment
%   same as e[n]
% tvLPCCoefficients: matrix of tvlpc coefficeints a[n,p]
%   LPC coefficient polynomial sampled at sampling instants 0:1:N-1

function [tvLPCParameters tvLPCResidual tvLPCCoefficients]=tvLPC(currentFrame,parameters,prevFrame,method)

%% handle inputs
    if(nargin<1)
        error('no inputs');
    elseif(nargin<2)
        params=[];
        prevFrame=zeros(windowSize,1);
        method='LS';
        disp('using default parameters and zero memory');
    elseif(nargin<3)
        prevFrame=zeros(windowSize,1);
        method='LS';
        disp('using zero memory');
    elseif(nargin<4)
        method='LS';
        disp('using least squares');
    end
    
    %% set default values
    if(~isfield(parameters,'samplingRate'))
        parameters.samplingRate=8000;
    end
    if(~isfield(parameters,'nPoles'))
        parameters.nPoles=round(parameters.samplingRate/1000)+2;
    end
    if(~isfield(parameters,'ord'))
        parameters.ord=3;
    end
    
%% Compute LPC coefficients   
    nPoles=parameters.nPoles;                                       % Number of poles
    nOrd=parameters.ord;                                            % Polynomial order
    windowSize=length(currentFrame);                                % Size of the signal frame
    S=zeros(windowSize,nPoles*(nOrd+1));                            % The matrix of previous samples
    s=currentFrame(:);                                              % Signal Frame
        
    for k=0:nOrd
        basis=([0:1:windowSize-1]/windowSize)'.^k;                  % n^k for 0<=n<=N-1
        for p=1:nPoles
            temp=[prevFrame(end-p+1:end);currentFrame(1:end-p)];    % s[n-p].*n^k
            S(:,k*nPoles+p)=temp.*basis;                            
        end        
    end    
    
    if(strcmp(method,'LS'))
        tvLPCParameters=S\s;                                            % Solve Least Squares problem
    else
        W=eye(windowSize);
        epsilon=1e-3;
        for i=1:3
            cvx_begin quiet
                variable a(nPoles*(nOrd+1),1)
                variable r(windowSize,1)
                minimize norm(W*r,1)
                subject to 
                r==s-S*a
            cvx_end
            residual=s-S*a;
            W=diag(1./(abs(residual)+epsilon));
        end        
%         cvx_begin quiet
%             variable a(nPoles*(nOrd+1),1)
%             variable r(windowSize,1)
%             minimize norm(r,1)
%             subject to 
%             r==s-S*a
%         cvx_end
        tvLPCParameters=a;clear a;
    end        
    tvLPCResidual=s-S*tvLPCParameters;
    
%% sample the LPC coefficient polynomial

    tvLPCCoefficients=zeros(windowSize,nPoles);        
    for p=1:nPoles
        temp=zeros(windowSize,1);
        for k=0:nOrd
            basis=([0:1:windowSize-1]/windowSize)'.^k;              % n^k
            temp=temp+tvLPCParameters(k*nPoles+p)*basis;            % { sum( a_pk * n^k ) }_k
        end
        tvLPCCoefficients(:,p)=temp;                                
    end        
    
end
