addpath('Routines')
addpath('Algorithms')

% -----------------------------------------------------------
%   Description
% -----------------------------------------------------------

% Monte Carlo simulation.
% Plots the probability of sub-meter accuracy VS. gain and phase mismatchs in all antennas.
% Location of the source, noise and the multipath channel are all random.


%% ----------------------------------------------------------
%   Parameters scenario
% -----------------------------------------------------------

% Number of antennas at all base stations
Nantennas = 70*ones(4,1);

% Carrier frequency [Hz]
carrierFreq = 7e9;

% Bandwidth [Hz]
bandwidth = 30e6;

% Width of the map where the source and stations are located
mapWidth = 100;

% Locations of the stations
baseStations = [45,45; 45,-45; -45,45; -45,-45];

% Multipath parameters
clusterDecayTime = 34e-9;
rayDecayTime = 29e-9;
angleStdDev = 26/180*pi;
clusterInterTime = 17e-9;
rayInterTime = 5e-9;

clusterDecayTime = 34e-9;
rayDecayTime = 29e-9;
angleStdDev = 26/180*pi;
clusterInterTime = 17e-6;
rayInterTime = 5e-6;

% ----------------------------------------------------------
%   Generic Parameters
% -----------------------------------------------------------

% Ration between LOS pulse energy and noise spectral density at each antennas [dB] (select multiple values)
EdivN0 = 10;

% Calibration errors
phaseStdDev = 0:10:100; % [�]
gainStdDevDB = 0:1.5:15; % [dB]

% Number of MC runs
MonteCarloRuns = 1e2;

% Oversampling factor
oversamplingFactor = 3;

% -----------------------------------------------------------
%   Parameters Algorithms
% -----------------------------------------------------------

% Starting resolution for grid based techniques [m].
startCellSize = 5;


%% ----------------------------------------------------------
%    Start simulation
% -----------------------------------------------------------

Nstations = length(Nantennas);
errorsPos = nan(length(phaseStdDev),length(gainStdDevDB),MonteCarloRuns);

hWait = waitbar(0,'Take a break, have a coffee, it make take some time...');
for p=1:MonteCarloRuns
    
    errorsPosRun = nan(length(phaseStdDev),length(gainStdDevDB));
    
    % Generate uniform circular arrays
    antLocs = cell(Nstations,1);
    fraunhoferDist = zeros(Nstations,1);
    for i=1:Nstations
        angleSep = 2*pi/Nantennas(i);
        radiusArray = 3e8/carrierFreq/4/sin(pi/Nantennas(i));
        orienation = rand*2*pi;
        antLocs{i} = radiusArray*[cos(angleSep*(1:Nantennas(i)).'+orienation),sin(angleSep*(1:Nantennas(i)).'+orienation)];
        antLocs{i} = ones(Nantennas(i),1)*baseStations(i,:) +antLocs{i};
        fraunhoferDist(i) = 2*(2*radiusArray)^2/3e8*carrierFreq;
    end

    % Location of the user
    mobileUser = .5*mapWidth*(rand(1,2)-.5);
    while any(sqrt(sum((baseStations-ones(Nstations,1)*mobileUser).^2,2))<fraunhoferDist)
        mobileUser = .5*mapWidth*(rand(1,2)-.5);
    end
    
    for s=1:length(phaseStdDev)
        for i=1:length(gainStdDevDB)
        
            % Computer powers
            receivedPowers = ones(Nstations,1);
            noisePower = 1/10^(EdivN0/10);

            % Generate multipath channels
            [trueTOAs, trueAOAs, trueAmplis ] = generateChannels( mobileUser, antLocs, receivedPowers, clusterDecayTime, rayDecayTime, clusterInterTime, rayInterTime, angleStdDev );

            % Apply MF on all antennas
            Nsamples = ceil(sqrt(2)*mapWidth/3e8*bandwidth*oversamplingFactor +20*clusterDecayTime*bandwidth*oversamplingFactor +2*oversamplingFactor); % Observations are taken until the energy of the received signal is negligible.
            [receivedSignals, receivedSignalsAfterMF, startIndex ] = generateUncalibratedSignals( antLocs, carrierFreq, bandwidth, oversamplingFactor, noisePower, Nsamples, trueTOAs, trueAOAs, trueAmplis, phaseStdDev(s), gainStdDevDB(i) );

            % Estimate TOA through threshold-based MF (positively biased)
            [TOAbyThresMF,~] = globalThresMF(receivedSignalsAfterMF,bandwidth,oversamplingFactor,startIndex,noisePower);

            % Estimate TOA's by finding the time a threshold is first corssed
            [TOAbyThres1stCross,snapshots] = globalFirstCrossThresMF(receivedSignalsAfterMF,bandwidth,oversamplingFactor,startIndex,noisePower);

            % Discard data from BS's where no signal is detected
            select = ~cellfun(@isempty,snapshots);
            if sum(select)>0
                snapshots = snapshots(select);
                TOAbyThresMF = TOAbyThresMF(select);
                TOAbyThres1stCross = TOAbyThres1stCross(select);
                activeBSLocs = baseStations(select,:);
                activeSensorsLocs = antLocs(select);
                
                % Estimate AOA's
                % cellfun(@(x)x(1)/pi*180,trueAOAs)
                AOAs = beamforming( activeSensorsLocs, snapshots, carrierFreq);

                mobileUserLocEstim = DiSouL( activeSensorsLocs, snapshots, noisePower, startCellSize, carrierFreq, bandwidth, TOAbyThresMF, mapWidth );
%                 mobileUserLocEstim = DPD( antLocs, receivedSignals, carrierFreq, bandwidth, oversamplingFactor, mapWidth, startCellSize, 5 );
%                 mobileUserLocEstim = Stansfield( activeBSLocs, TOAbyThres1stCross, AOAs );

                if ~isempty(mobileUserLocEstim) && ~any(isnan(mobileUserLocEstim))
                    errorsPosRun(s,i) = norm(mobileUserLocEstim-mobileUser);
                end
            end
        end
    end
    errorsPos(:,:,p) = errorsPosRun;    
    
    waitbar(p/MonteCarloRuns,hWait);
end
close(hWait);


%% ----------------------------------------------------------
%   Output Results
% -----------------------------------------------------------

figure
probab = mean(errorsPos<=1,3);
[X,Y] = meshgrid(gainStdDevDB,phaseStdDev);
meshc(X,Y,probab);
ylabel('Phase standard deviation [�]');
xlabel('Gain standard deviation [dB]');

figure
fitobject = fit([X(:),Y(:)],probab(:),'loess');
plot(fitobject,[X(:),Y(:)],probab(:))
smoothedZ = fitobject(X(:),Y(:));

%contourf(gainStdDevDB,phaseStdDev,probab);
%contourc(gainStdDevDB,phaseStdDev,probab);

print('../results/Fig9','-dpng')