addpath('Routines')
addpath('Algorithms')

% -----------------------------------------------------------
%   Description
% -----------------------------------------------------------

% Monte Carlo simulation.
% Plots root mean squared error (rMSE) and probability of sub-meter accuracy VS. E/N0.
% Location of the source, noise and the multipath channel are all random.


%% ----------------------------------------------------------
%   Parameters scenario
% -----------------------------------------------------------

% Number of antennas at all base stations
numAntennas = 70*ones(4,1);

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

% Ration between LOS pulse energy and noise spectral density at each antenna [dB] (select multiple values)
EdivN0 = -15:5:25;

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

Nstations = length(numAntennas);
errorsPos = nan(length(algorithm),length(EdivN0),MonteCarloRuns);

hWait = waitbar(0,'Take a break, have a coffee, it make take some time...');
for p=1:MonteCarloRuns
    
    errorsPosRun = nan(length(algorithm),length(EdivN0));
    
    % Generate uniform circular arrays
    antLocs = cell(Nstations,1);
    fraunhoferDist = zeros(Nstations,1);
    for i=1:Nstations
        angleSep = 2*pi/numAntennas(i);
        radiusArray = 3e8/carrierFreq/4/sin(pi/numAntennas(i));
        orienation = rand*2*pi;
        antLocs{i} = radiusArray*[cos(angleSep*(1:numAntennas(i)).'+orienation),sin(angleSep*(1:numAntennas(i)).'+orienation)];
        antLocs{i} = ones(numAntennas(i),1)*baseStations(i,:) +antLocs{i};
        fraunhoferDist(i) = 2*(2*radiusArray)^2/3e8*carrierFreq;
    end
    
    % Location of the source
    mobileUser = .5*mapWidth*(rand(1,2)-.5);
    while any(sqrt(sum((baseStations-ones(Nstations,1)*mobileUser).^2,2))<fraunhoferDist)
        mobileUser = .5*mapWidth*(rand(1,2)-.5);
    end
    
    for s=1:length(EdivN0)
        
        % Computer powers
        receivedPowers = ones(Nstations,1);
        noisePower = 1/10^(EdivN0(s)/10);
        
        % Generate multipath channels
        [trueTOAs, trueAOAs, trueAmplis ] = generateChannels( mobileUser, antLocs, receivedPowers, clusterDecayTime, rayDecayTime, clusterInterTime, rayInterTime, angleStdDev );

        % Apply MF on all antennas
        Nsamples = ceil(sqrt(2)*mapWidth/3e8*bandwidth*oversamplingFactor +20*clusterDecayTime*bandwidth*oversamplingFactor +2*oversamplingFactor); % Observations are taken until the energy of the received signal is negligible.
        [receivedSignals, receivedSignalsAfterMF, startIndex ] = generateSignals( antLocs, carrierFreq, bandwidth, oversamplingFactor, noisePower, Nsamples, trueTOAs, trueAOAs, trueAmplis );
        
        % Estimate TOA through threshold-based MF (positively biased)
        [TOAbyThresMF,~] = globalThresMF(receivedSignalsAfterMF,bandwidth,oversamplingFactor,startIndex,noisePower);

        % Estimate TOA's by finding the time a threshold is first corssed
        [TOAbyThres1stCross,snapshots] = globalFirstCrossThresMF(receivedSignalsAfterMF,bandwidth,oversamplingFactor,startIndex,noisePower);
        
        % Discard data from BS's where no signal is detected
        select = ~cellfun(@isempty,snapshots);
        if sum(select)>0
            snapshots = snapshots(select);
            activeTOAbyThresMF = TOAbyThresMF(select);
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
                        mobileUserLocEstim = SR_LS( activeBSLocs, activeTOAbyThresMF );
                    case 31
                        mobileUserLocEstim = Stansfield( activeBSLocs, activeTOAbyThresMF, AOAs );
                    case 41 % DPD
                        mobileUserLocEstim = DPD( antLocs, receivedSignals, carrierFreq, bandwidth, oversamplingFactor, mapWidth, startCellSize, 5 );
                    case 42 % DiSouL
                        [mobileUserLocEstim,losBS] = DiSouL( activeSensorsLocs, snapshots, noisePower, startCellSize, carrierFreq, bandwidth, activeTOAbyThresMF, mapWidth );
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
        case 41
            textLegend = [textLegend, '''DPD'''];
        case 42
            textLegend = [textLegend, '''DiSouL'''];
    end
    if a~=length(algorithm)
        textLegend = [textLegend, ','];
    end
end 

markerTypes = '<o+s^hxd>p*.v';
colorTypes = 'bgrymc';

%rMSE = sqrt(mean(errorsPos.^2,3));
%semilogy(EdivN0,rMSE.');
%xlabel('SNR per antenna [dB]');
%ylabel('rMSE [m]');
%eval(['legend(',textLegend,');']);
%ylim([0,mapWidth]);

% hold
% semilogy(EdivN0,sqrt(mean(CRB,2)));

figure
probab = sum(errorsPos<1,3)/MonteCarloRuns;
plot(EdivN0,probab.');
xlabel('SNR [dB]');
ylabel('Probability of submeter accuracy');
ylim([0,1]);
eval(['legend(',textLegend,');']);
print('../results/Fig5','-dpng')