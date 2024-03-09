addpath('Routines')
addpath('Algorithms')

% -----------------------------------------------------------
%   Description
% -----------------------------------------------------------

% Monte Carlo simulation.
% Plots root mean squared error (rMSE) and probability of sub-meter accuracy.
% Location of the source, configuration of the arrays and the multipath channel is random.


%% ----------------------------------------------------------
%   Parameters scenario
% -----------------------------------------------------------

% Number of sources and sensors
Nantennas = 70*ones(4,1);

% Carrier frequency [Hz]
carrierFreq = 7e9;

% Bandwidth [Hz]
bandwidth = 30e6;

% Width of the map where emitters and sensors are located
mapWidth = 100;

% Locations of the emitters
mobileUser = .5*mapWidth*(rand(1,2)-.5);

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
%   - 100 Test algorithm

algorithm = [11,21,31,41,42];

% Ration between LOS pulse energy and noise spectral density at each antenna [dB] (select multiple values)
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
errorsPos = nan(length(algorithm),MonteCarloRuns);
elapsedTime = nan(length(algorithm),MonteCarloRuns);
adjAmplis = cell(MonteCarloRuns,1);

% Computer powers
receivedPowers = ones(Nstations,1);
noisePower = 1/10^(EdivN0/10);

hWait = waitbar(0,'Take a break, have a coffee, it make take some time...');
for p=1:MonteCarloRuns
    
    errorsPosRun = nan(length(algorithm),1);
    elapsedTimeRun = nan(length(algorithm),1);

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
    
    % Location of the emitter
    mobileUser = .5*mapWidth*(rand(1,2)-.5);
    while any(sqrt(sum((baseStations-ones(Nstations,1)*mobileUser).^2,2))<fraunhoferDist)
        mobileUser = .5*mapWidth*(rand(1,2)-.5);
    end

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

        % Find the signal strenghts in the snapshot
        adjAmplisRun = findSnapshotAmplis( bandwidth, trueTOAs(select), trueAmplis(select), TOAbyThres1stCross );
        adjAmplisRun = cellfun(@(x)abs(x).^2,adjAmplisRun,'UniformOutput',false);

        % Estimate AOA's
        AOAs = beamforming( activeSensorsLocs, snapshots, carrierFreq);

        % Run localization algorithms
        for a=1:length(algorithm)
            tic
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
                case 100 % Test algorithm
                    mobileUserLocEstim = test_alg( mobileUser, activeSensorsLocs, carrierFreq, bandwidth, noisePower, cellfun(@(x)x(1),trueTOAs(select)), cellfun(@(x)x(1),trueAOAs(select)), cellfun(@(x)x(1),trueAmplis(select)), TOAbyThres1stCross );
            end
            elapsedTimeRun(a) = toc;

            if ~isempty(mobileUserLocEstim) && ~any(isnan(mobileUserLocEstim))
                errorsPosRun(a) = norm(mobileUserLocEstim-mobileUser);
            end
        end

    end
        
    errorsPos(:,p) = errorsPosRun;
    elapsedTime(:,p) = elapsedTimeRun;
    adjAmplis{p} = adjAmplisRun;
    
    waitbar(p/MonteCarloRuns,hWait);
end
close(hWait);


%% ----------------------------------------------------------
%   Output Results
% -----------------------------------------------------------

adjAmplis = vertcat(adjAmplis{:});
numArrivals99 = zeros(length(adjAmplis),1);
numArrivals999 = zeros(length(adjAmplis),1);
numArrivalsPower = zeros(length(adjAmplis),1);
for i=1:length(adjAmplis)
    orderedPows = sort(adjAmplis{i},'descend');
    numArrivals99(i) = find(cumsum(orderedPows)>=.99*sum(orderedPows),1);
    numArrivals999(i) = find(cumsum(orderedPows)>=.999*sum(orderedPows),1);
    numArrivalsPower(i) = sum(adjAmplis{i}*mean(Nantennas)>noisePower);
end
disp(['The average number of arrivals containing 99% of the energy is ', num2str(mean(numArrivals99)), '.'])
disp(['The average number of arrivals containing 99.9% of the energy is ', num2str(mean(numArrivals999)), '.'])
disp(['The average number of arrivals whose power is larger than the aggregated noise power is ', num2str(mean(numArrivalsPower)), '.'])

textLegend = {};
avTime = mean(elapsedTime,2);
for a=1:length(algorithm)
    switch algorithm(a)
        case 11
            textToAdd = 'IV';
        case 21
            textToAdd = 'SR-LS';
        case 31
            textToAdd = 'Stansfield';
        case 41
            textToAdd = 'DPD';
        case 42
            textToAdd = 'DiSouL';
        case 100
            textToAdd = 'Test alg';
    end
    textLegend = [textLegend,{textToAdd}];
    disp(['The average execution time for ', textToAdd, ' is ', num2str(mean(elapsedTime(a,:))*1e3), ' ms.'])
end 

markerTypes = '<o+s^hxd>p*.v';
colorTypes = 'bgrymc';

%figure
%rMSE = sqrt(mean(errorsPos.^2,2));
%bar(rMSE);
%ylabel('rMSE [m]');
%set(gca,'XTickLabel',textLegend);

%figure
%probab = sum(errorsPos<1,2)/MonteCarloRuns;
%bar(probab);
%ylabel('Probability of submeter accuracy');
%ylim([0,1]);
%set(gca,'XTickLabel',textLegend);

figure
hold on
for a=1:length(algorithm)
    cdfplot(errorsPos(a,:))
end
hold off
xlim([0,mapWidth/2]);
legend(textLegend);

print('../results/Fig4','-dpng')