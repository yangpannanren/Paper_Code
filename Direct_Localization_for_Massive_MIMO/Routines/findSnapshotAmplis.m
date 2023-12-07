function [ adjCompAmplis ] = findSnapshotAmplis( bandwidth, TOAs, compAmplis, samplInstants )
%DIRECT Summary of this function goes here
%   Detailed explanation goes here

% Coded by Nil Garcia

%   Inputs:
%   - bandwidth [Hz]
%       The bandwidth of the emitted signal.
%
%   Outputs:
%   - adjCompAmplis
    

% -----------------------------------------------------------
%   Check that the inputs comply with the correct format.
% -----------------------------------------------------------
    
    % Check 'bandwidth'
    if ~isscalar(bandwidth)
        error('Localization error: The sampling frequency is not a scalar.')
    elseif ~isreal(bandwidth)
        error('Localization error: The sampling frequency is not a real number.')
    elseif ~isfinite(bandwidth)
        error('Localization error: Infinite is not valid for the sampling frequency.')
    elseif bandwidth<=0
        error('Localization error: The sampling frequency must be positive.')
    end
    
    numBS = length(TOAs);
    

% -----------------------------------------------------------
%   Find amplitudes of the arrivals in the snapshots
% -----------------------------------------------------------
    
%     adjCompAmplis = cellfun(@(x,y,z)exp(-((x-y)*bandwidth*pi).^2./4./log(2)).*z,TOAs,num2cell(samplInstants),compAmplis,'UniformOutput','false');

    adjCompAmplis = cell(numBS,1);
    for i=1:numBS
        adjCompAmplis{i} = exp(-((TOAs{i}-samplInstants(i))*bandwidth*pi).^2./4./log(2)).*compAmplis{i};
    end
 
end