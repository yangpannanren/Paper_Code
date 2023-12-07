function [ steeringVectors, diffSteeringVectors ] = angleToSteering( antennaLocs, angles, carrierFreq, varargin )

% Coded by Nil Garcia

%   Inputs:r
%   - antennaLocs [m]
%       Matrix of size (Nantennas)x2 matrix which stacks the
%       positions of the sensors.
%   - angles [rad]
%       Vector of length (numAngles)
%       Angles of the incidencent signals. The positive part of the x-axis
%       corresponds to 0 rad.
%   - carrierFreq [Hz]
%       The carrier frequency of the emitted signals.
%
%   Outputs:
%   - steeringVectors
%       Matrix of dimensions (Nantennas)x(numAngles)
%       Each entry corresponds to the phase offset for a sensor for a
%       signal arriving from a certain angle.

    Nantennas = size(antennaLocs,1);
    waveLength = 3e8/carrierFreq;
    centerOfGravity = mean(antennaLocs,1);
    
    steeringVectors = exp((antennaLocs-ones(Nantennas,1)*centerOfGravity)*[cos(angles(:)),sin(angles(:))].'*2*pi/waveLength*1i);
    diffSteeringVectors = 2*pi/waveLength*1i*((antennaLocs-ones(Nantennas,1)*centerOfGravity)*[-sin(angles(:)),cos(angles(:))].').*steeringVectors;
    
    if length(varargin)==2
        phaseStdDev = varargin{1};
        gainStdDevDB = varargin{2};
%         ampliStdDev = 10^(varargin{2}/20);
        
        phases = phaseStdDev*randn(Nantennas,1)/180*pi;
%         gains = lognrnd(-ampliStdDev^2/2,ampliStdDev,[Nantennas,1]);
        gainsDB = -log(10)/40*gainStdDevDB^2 +gainStdDevDB*randn(Nantennas,1);
        gains = 10.^(gainsDB./20);
        
        steeringVectors = diag(gains.*exp(1i*phases))*steeringVectors;
        diffSteeringVectors = diag(gains.*exp(1i*phases))*diffSteeringVectors;
    end


end
