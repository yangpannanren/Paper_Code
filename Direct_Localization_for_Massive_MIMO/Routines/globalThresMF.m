function [ TOAs, snapshot ] = globalThresMF(receivedSignals,bandwidth,oversamplingFactor,startIndex,noisePower)
%TOABYMG Summary of this function goes here
%   Estimates the TOA at multiple arrays based on simple match filter plus
%   largest peak.
    
% Coded by Nil Garcia

%   Inputs:
%   - receivedSignal
%       Cell array of dimensions (numBS)x1
%       The m-th cell contains the (numSensors(m))x1 samples from the
%       single snapshot.
%   - bandwidth [Hz]
%       The sampling frequency at the sensors.
%   - oversamplingFactor
%       The ratio between the sampling frequency and the bandwidth.
%   - startIndex
%       Sample index of 'receivedSignal' corresponding to delay 0.
%   - noisePower
%
%   Outputs:
%   - TOAs [s]
%       Vector of length (numBs) with the estimated TOAs.
%   - snapshot

    % Check 'receivedSignals'
    if ~iscell(receivedSignals)
        error('TOA estimation error: The received signals must be represented by cell array.')
    elseif ~isvector(receivedSignals)
        error('TOA estimation error: The received signals must be a vector cell array.')
    elseif any(~cellfun(@ismatrix,receivedSignals))
        error('TOA estimation error: The received signals for each array must be represented by a matrix.')
    elseif any(~cellfun(@isnumeric,receivedSignals))
        error('TOA estimation error: The received signals coordinates must be numbers.')
    end
    numBS = length(receivedSignals);
    numSensors = cellfun(@(x)size(x,1),receivedSignals);
    numSamples = size(receivedSignals{1},2);
    if any(cellfun(@(x)size(x,2),receivedSignals)~=numSamples)
        error('TOA estimation error: The received signals at all base stations must be of equal length.')
    end
    
    % Check 'bandwidth'
    if ~isscalar(bandwidth)
        error('TOA estimation error: The bandwidth is not a scalar.')
    elseif ~isreal(bandwidth)
        error('TOA estimation error: The bandwidth is not a real number.')
    elseif ~isfinite(bandwidth)
        error('TOA estimation error: Infinite is not a valid value for the bandwidth.')
    elseif bandwidth<=0
        error('TOA estimation error: The bandwidth must be positive.')
    end
    
    % Check 'startIndex'
    if ~isscalar(startIndex)
        error('TOA estimation error: The starting index is not a scalar.')
    elseif ~isreal(startIndex) || startIndex-round(startIndex)~=0
        error('TOA estimation error: The starting index is not an integer.')
    elseif ~isfinite(startIndex)
        error('TOA estimation error: Infinite is not valid for the starting index.')
    elseif startIndex<=0
        error('TOA estimation error: The starting index must be positive.')
    end
    
    % Check 'oversamplingFactor'
    if ~isscalar(oversamplingFactor)
        error('TOA estimation error: The oversampling factor is not a scalar.')
    elseif ~isreal(oversamplingFactor) || oversamplingFactor-round(oversamplingFactor)~=0
        error('TOA estimation error: The oversampling factor is not an integer.')
    elseif ~isfinite(oversamplingFactor)
        error('TOA estimation error: Infinite is not valid for the oversampling factor.')
    elseif oversamplingFactor<=0
        error('TOA estimation error: The oversampling factor must be positive.')
    end
    
    % Check 'noisePower'
    if ~isscalar(noisePower)
        error('TOA estimation error: The noise power is not a scalar.')
    elseif ~isreal(noisePower)
        error('TOA estimation error: The noise power is not a real number.')
    elseif ~isfinite(noisePower)
        error('TOA estimation error: Infinite is not a valid value for the noise power.')
    elseif noisePower<0
        error('TOA estimation error: The noise power cannot be negative.')
    end


% -----------------------------------------------------------
%   Compute threshold
% -----------------------------------------------------------
    
    thres = zeros(numBS,1);
    if noisePower>0
        for l=1:numBS
            numSlots = (numSamples-startIndex+1)/oversamplingFactor;    
            candidateThres = 0:numSensors(l)*noisePower/10:20*numSensors(l)*noisePower;
            probFA = 1e-2;
            
            qNoise = chi2cdf(candidateThres/noisePower*2,2*numSensors(l),'upper');
            temp = candidateThres(1+((1-qNoise).^numSlots-1)./numSlots./qNoise < probFA);
            thres(l) = temp(1);
        end
    end
    
    
% -----------------------------------------------------------
%   Estimate TOA's
% -----------------------------------------------------------

    TOAs = nan(numBS,1);
    snapshot = cell(numBS,1);
    warning('off','signal:findpeaks:largeMinPeakHeight')
    for l=1:numBS
        
        NCsignal = sum(abs(receivedSignals{l}).^2,1);
        [~,ind] = findpeaks(NCsignal,'MINPEAKHEIGHT',thres(l),'SORTSTR','none');
        
        if ~isempty(ind)
            % Find peak location by parabolic-fit (see "Methods for estimation of subsample time delays of digitized echo signals" for more details)
            peakOffset = (sqrt(NCsignal(ind(1)-1))-sqrt(NCsignal(ind(1)+1)))/2/(sqrt(NCsignal(ind(1)-1))-2*sqrt(NCsignal(ind(1)))+sqrt(NCsignal(ind(1)+1)));
            peakInd = ind(1)+peakOffset;
            
            TOAs(l) = (peakInd-startIndex)/oversamplingFactor/bandwidth;
            snapshot{l} = receivedSignals{l}(:,ind(1));
        end
    end
    warning('on','signal:findpeaks:largeMinPeakHeight')
    
    
end
