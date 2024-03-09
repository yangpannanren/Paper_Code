function [ receivedSignals, signalAfterMF, startIndex ] = generateUncalibratedSignals( sensorsLocs, carrierFreq, bandwidth, oversamplingFactor, noisePower, numSamples, TOAs, AOAs, compAmpli, phaseStdDev, gainStdDevDB )

    
% The emitted signals are narrow-band Gaussian pulses.

% Coded by Nil Garcia

%   Inputs:
%   - sensorsLocs [m]
%       Cell vector array of length (numBS)
%       The m-th cell contains l (numSensors(m))x2 matrix which stacks the
%       positions of the sensors within such array.
%   - carrierFreq [Hz]
%       The carrier frequency of the emitted signal.
%   - bandwidth [Hz]
%       The bandwidth of the emitted signal.
%   - oversamplingFactor
%       The oversampling factor at reception
%   - noisePower
%       The noise power
%   - numSamples
%       Number of observations
%   - TOAs [s]
%       Cell vector array of length (numBS)
%       The m-th cell contains a vector of all TOAs at the m-th BS.
%   - AOAs [radians]
%       Cell vector array of length (numBS)
%       The m-th cell contains a vector of all AOAs at the m-th BS.
%   - compAmpli
%       Cell vector array of length (numBS)
%       The m-th cell contains a vector of all complex amplitudes at the m-th BS.
%
%   Outputs:
%   - signalAfterMF
%       Cell array of dimensions (numBS)x1
%       The m-th cell contains the (numSensors(m))x(numSamples) received signals after MF.
%   - startIndex
%       Sample index of 'receivedSignals' corresponding to delay 0.
    

% -----------------------------------------------------------
%   Check that the inputs comply with the correct format.
% -----------------------------------------------------------
    
    % Check 'sensorsLocs'
    if ~iscell(sensorsLocs)
        error('Signal generation error: The sensors locations must be represented by cell array.')
    elseif ~isvector(sensorsLocs)
        error('Signal generation error: The sensors locations must be a vector cell array.')
    elseif any(~cellfun(@ismatrix,sensorsLocs))
        error('Signal generation error: The sensors locations for each array must be represented by a matrix.')
    elseif any(cellfun(@(x)size(x,2),sensorsLocs)~=2)
        error('Signal generation error: The number of coordinates of the sensors locations must be 2.')
    elseif any(~cellfun(@isreal,sensorsLocs))
        error('Signal generation error: The sensors locations coordinates must be real numbers.')
    end
    numBS = length(sensorsLocs);
    numSensors = cellfun(@(x)size(x,1),sensorsLocs);

    % Check 'carrierFreq'
    if ~isscalar(carrierFreq)
        error('Signal generation error: The carrier frequeny is not a scalar.')
    elseif ~isreal(carrierFreq)
        error('Signal generation error: The carrier frequeny is not a real number.')
    elseif ~isfinite(carrierFreq)
        error('Signal generation error: Infinite is not a valid value for the carrier frequeny.')
    elseif carrierFreq<=0
        error('Signal generation error: The carrier frequeny must be positive.')
    end
    
    % Check 'bandwidth'
    if ~isscalar(bandwidth)
        error('Signal generation error: The bandwidth is not a scalar.')
    elseif ~isreal(bandwidth)
        error('Signal generation error: The bandwidth is not a real number.')
    elseif ~isfinite(bandwidth)
        error('Signal generation error: Infinite is not a valid value for the bandwidth.')
    elseif bandwidth<=0
        error('Signal generation error: The bandwidth must be positive.')
    end
    
    % Check 'oversamplingFactor'
    if ~isscalar(oversamplingFactor)
        error('Signal generation error: The oversampling factor is not a scalar.')
    elseif ~isreal(oversamplingFactor) || oversamplingFactor-round(oversamplingFactor)~=0
        error('Signal generation error: The oversampling factor is not an integer.')
    elseif ~isfinite(oversamplingFactor)
        error('Signal generation error: Infinite is not valid for the oversampling factor.')
    elseif oversamplingFactor<=0
        error('Signal generation error: The oversampling factor must be positive.')
    end
    
    % Check 'numSamples'
    if ~isscalar(numSamples)
        error('Signal generation error: The number of samples factor is not a scalar.')
    elseif ~isreal(numSamples) || numSamples-round(numSamples)~=0
        error('Signal generation error: The number of samples is not an integer.')
    elseif ~isfinite(numSamples)
        error('Signal generation error: Infinite is not valid for the number of samples.')
    elseif numSamples<=0
        error('Signal generation error: The number of samples must be positive.')
    end
    
    % Check 'noisePower'
    if ~isscalar(noisePower)
        error('Signal generation error: The noise power is not a scalar.')
    elseif ~isreal(noisePower)
        error('Signal generation error: The noise power is not a real number.')
    elseif ~isfinite(noisePower)
        error('Signal generation error: Infinite is not a valid value for the noise power.')
    elseif noisePower<0
        error('Signal generation error: The noise power cannot be negative.')
    end
    
    
% -----------------------------------------------------------
%   Compute the received signals and apply match filter
% -----------------------------------------------------------
    
    receivedSignals = cell(numBS,1);
    
    offsetSamples = ceil(3*oversamplingFactor/pi*sqrt(2*log(2)));
    
    for l=1:numBS

        receivedSignals{l} = sqrt(noisePower/2)*randn(numSensors(l),numSamples) +1i*sqrt(noisePower/2)*randn(numSensors(l),numSamples);

        % Simulate signal at receivers
        for i=1:length(TOAs{l})
            receivedSignals{l} = receivedSignals{l} +compAmpli{l}(i)...
                *angleToSteering(sensorsLocs{l},AOAs{l}(i),carrierFreq,phaseStdDev,gainStdDevDB)...
                *oversamplingFactor^(-1/2)*(pi/log(2))^(1/4)...
                *exp(-((0:numSamples-1)./oversamplingFactor-bandwidth*TOAs{l}(i)-offsetSamples/oversamplingFactor).^2.*pi.^2./2./log(2));
        end
    end
    
    waveform = exp(-((-offsetSamples:offsetSamples)./oversamplingFactor.*pi).^2./2./log(2));
    waveform = waveform./norm(waveform);

    signalAfterMF = cell(numBS,1);
    for l=1:numBS
        signalAfterMF{l} = zeros(numSensors(l),numSamples+length(waveform)-1);
        for s=1:numSensors(l)
            signalAfterMF{l}(s,:) = conv(receivedSignals{l}(s,:),waveform,'full');
        end
    end
    
    % Index of the sample corresponding to zero delay.
    startIndex = length(waveform);
    
%     % TEST
%     plot((0:length(signalAfterMF{1}(1,startIndex:end))-1)/bandwidth/oversamplingFactor*1e6,abs(signalAfterMF{1}(1,startIndex:end)))
%     xlabel('us')
%     % TEST
end
