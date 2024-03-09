addpath('Routines')
addpath('Algorithms')

% -----------------------------------------------------------
%   Description
% -----------------------------------------------------------

% Monte Carlo simulation.
% Plots root probability of sub-meter accuracy VS. bandwidth of the signal.
% Location of the source, noise and the multipath channel are all random.


%% ----------------------------------------------------------
%   Parameters scenario
% -----------------------------------------------------------

% Number of sources and sensors
Nantennas = 70*ones(4,1);

% Carrier frequency [Hz]
carrierFreq = 7e9;

% Bandwidth [Hz]
bandwidth = 10e6:10e6:100e6;

% Width of the map where emitters and sensors are located
mapWidth = 100;

% Locations of the arrays of sensors
baseStations = [45,45; 45,-45; -45,45; -45,-45];

% Multipath parameters
clusterDecayTime = 34e-9;
rayDecayTime = 29e-9;
clusterInterTime = 17e-9;
rayInterTime = 5e-9;
angleStdDev = 26/180*pi;


% ----------------------------------------------------------
%   Generic Parameters
% -----------------------------------------------------------

% Localization technique selection (only one)
%   - 11 Instrumental Variables estimator (IV)
%   - 21 Squared-range-based Least Squares (SR-LS)
%   - 31 Stansfield estimator
%   - 41 Direct Position Determination (DPD)
%   - 42 DiSouL

algorithm = [11,21,31,41,42];

% Ration between LOS pulse energy and noise spectral density at each antennas [dB] (select multiple values)
EdivN0 = 10;

% Number of particles per E/N0 value
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
errorsPos = nan(length(algorithm),length(bandwidth),MonteCarloRuns);

hWait = waitbar(0,'Take a break, have a coffee, it make take some time...');
for p=1:MonteCarloRuns
    
    errorsPosRun = nan(length(algorithm),length(bandwidth));

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
    
    % Location of the source
    mobileUser = .5*mapWidth*(rand(1,2)-.5);
    while any(sqrt(sum((baseStations-ones(Nstations,1)*mobileUser).^2,2))<fraunhoferDist)
        mobileUser = .5*mapWidth*(rand(1,2)-.5);
    end
    
    for s=1:length(bandwidth)
        
        % Computer powers
        receivedPowers = ones(Nstations,1);
        noisePower = 1/10^(EdivN0/10);

        % Generate multipath channels
        [trueTOAs, trueAOAs, trueAmplis ] = generateChannels( mobileUser, antLocs, receivedPowers, clusterDecayTime, rayDecayTime, clusterInterTime, rayInterTime, angleStdDev );

        % Apply MF on all antennas
        Nsamples = ceil(sqrt(2)*mapWidth/3e8*bandwidth(s)*oversamplingFactor +20*clusterDecayTime*bandwidth(s)*oversamplingFactor +2*oversamplingFactor); % Observations are taken until the energy of the received signal is negligible.
        [receivedSignals, receivedSignalsAfterMF, startIndex ] = generateSignals( antLocs, carrierFreq, bandwidth(s), oversamplingFactor, noisePower, Nsamples, trueTOAs, trueAOAs, trueAmplis );

        % Estimate TOA through threshold-based MF (positively biased)
        [TOAbyThresMF,~] = globalThresMF(receivedSignalsAfterMF,bandwidth(s),oversamplingFactor,startIndex,noisePower);

        % Estimate TOA's by finding the time a threshold is first corssed
        [TOAbyThres1stCross,snapshots] = globalFirstCrossThresMF(receivedSignalsAfterMF,bandwidth(s),oversamplingFactor,startIndex,noisePower);
        
        % Discard data from BS's where no signal is detected
        select = ~cellfun(@isempty,snapshots);
        if sum(select)>0
            snapshots = snapshots(select);
            TOAbyThresMF = TOAbyThresMF(select);
            TOAbyThres1stCross = TOAbyThres1stCross(select);
            activeBSLocs = baseStations(select,:);
            activeSensorsLocs = antLocs(select);

            % Estimate AOA's
            AOAs = beamforming( activeSensorsLocs, snapshots, carrierFreq);

            % Run localization algorithms
            for a=1:length(algorithm)
                switch algorithm(a)
                    case 11
                        mobileUserLocEstim = IV( activeBSLocs, AOAs );
                    case 21
                        mobileUserLocEstim = SR_LS( activeBSLocs, TOAbyThresMF );
                    case 31
                        mobileUserLocEstim = Stansfield( activeBSLocs, TOAbyThresMF, AOAs );
                    case 41 % DPD
                        mobileUserLocEstim = DPD( antLocs, receivedSignals, carrierFreq, bandwidth(s), oversamplingFactor, mapWidth, startCellSize, 5 );
                    case 42 % DiSouL
                        [mobileUserLocEstim,losBS] = DiSouL( activeSensorsLocs, snapshots, noisePower, startCellSize, carrierFreq, bandwidth(s), TOAbyThresMF, mapWidth );
                end

                if ~isempty(mobileUserLocEstim) && ~any(isnan(mobileUserLocEstim))
                    errorsPosRun(a,s) = norm(mobileUserLocEstim-mobileUser);
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

textLegend = [];
for a=1:length(algorithm)
    switch algorithm(a)
        case 11
            textLegend = [textLegend, '''IV'''];
        case 21
            textLegend = [textLegend, '''SR-LS'''];
        case 31
            textLegend = [textLegend, '''Stansfield'''];
        case 41 % DPD
            textLegend = [textLegend, '''DPD'''];
        case 42
            textLegend = [textLegend, '''DiSouL'''];
    end
    if a~=length(algorithm)
        textLegend = [textLegend, ','];
    end
end 

figure
probab = sum(errorsPos<1,3)/MonteCarloRuns;
plot(bandwidth/1e6,probab.');
xlabel('Bandwidth [MHz]');
ylabel('Probability of submeter accuracy');
ylim([0,1]);
eval(['legend(',textLegend,');']);

print('../results/Fig6','-dpng')
