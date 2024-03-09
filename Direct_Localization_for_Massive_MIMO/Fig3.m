addpath('Routines')
addpath('Algorithms')

% -----------------------------------------------------------
%   Description
% -----------------------------------------------------------

% Monte Carlo simulation.
% Plots probability of sub-meter accuracy VS. the SNR of the snapshots.
% Location of the source, configuration of the arrays and multipath components' strengths and phases are random.


%% ----------------------------------------------------------
%   Parameters scenario
% -----------------------------------------------------------

% Number of antennas at all base stations
Nantennas = 70*ones(4,1);

% Carrier frequency [Hz]
carrierFreq = 7e9;

% Width of the map where the source and stations are located
mapWidth = 100;

% Locations of the source
mobileUser = [18,31];

% Locations of the stations
baseStations = [45,45; 45,-45; -45,45; -45,-45];

% Position of the reflector
reflector = [25,-7];

% BSs receiving a NLOS path
BSnlos = [1,2,4];


% ----------------------------------------------------------
%   Generic Parameters
% -----------------------------------------------------------

% Square of the weight parameter
Nstations = length(Nantennas);
wSquare = .6:.2:5;

% SNR [dB] (select multiple values)
SNR = 0:10:20;

% Number of particles per E/N0 value
MonteCarloRuns = 1e2;


% -----------------------------------------------------------
%   Parameters Algorithms
% -----------------------------------------------------------

% Starting resolution for grid based techniques [m].
startCellSize = 5;


%% ----------------------------------------------------------
%    Start simulation
% -----------------------------------------------------------

errorsPos = nan(length(wSquare),length(SNR),MonteCarloRuns);

hWait = waitbar(0,'Take a break, have a coffee, it make take some time...');
for p=1:MonteCarloRuns
    
    errorsPosRun = nan(length(wSquare),length(SNR));

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
    
    for s=1:length(SNR)
        
        % Computer powers
        receivedPowers = 1;
        noisePower = 1/10^(SNR(s)/10);
            
        % Computing the received signals
        [snapshots,trueAOAs] = generateSignalsSimpleScenario( mobileUser, antLocs, receivedPowers, carrierFreq, mapWidth, reflector, BSnlos);
        
        % Add noise
        snapshots = cellfun(@(x,y)x+sqrt(noisePower/2)*(randn(y,1) +1i*randn(y,1)),snapshots,num2cell(Nantennas),'UniformOutput',false);

        for w=1:length(wSquare)

            % Run localization algorithms
            mobileUserEstim = DiSouLfixedWnoTOA( antLocs, snapshots, noisePower, startCellSize, carrierFreq, mapWidth, wSquare(w) );

            if ~isempty(mobileUserEstim) && ~any(isnan(mobileUserEstim))
                errorsPosRun(w,s) = norm(mobileUserEstim-mobileUser);
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

% Prob <1m error
probab = mean(errorsPos<1,3);
h = figure;
hold on
handlers = [];
for s=1:length(SNR)
    handlers = [handlers, plot(wSquare,probab(:,s),'DisplayName',['SNR = ',num2str(SNR(s)),' dB.'])];
end
hold off
ylabel('Probability of submeter accuracy');
xlabel('w^2')
ylim([0,1]);
legend(handlers);

print('../results/Fig3','-dpng')