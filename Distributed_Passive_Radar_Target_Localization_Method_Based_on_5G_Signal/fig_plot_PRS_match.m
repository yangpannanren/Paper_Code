function fig_plot_PRS_match()
% 函数功能：画5G_PRS匹配图
% 参考链接：https://ww2.mathworks.cn/help/5g/ug/nr-prs-positioning.html?searchHighlight=validateCarriers&s_tid=srchtitle_support_results_1_validateCarriers
rng('default');
nFrames = 1; % Number of 10 ms frames
fc = 5.1e9;  % Carrier frequency (Hz)
UEPos = [2e3 4e3 0]; % In meters. As per TR 38.901 Table 7.4.1-1,
                     % UE height for 'UMa' path loss scenario is >= 1.5 meters and <= 22.5 meters
numgNBs = 4;
gNBPos = {[0,4e3,0],[0,2e3,0],[2e3,2e3,0],[0,0,0]};
cellIds = randperm(1008,numgNBs) - 1;
% Configure carrier properties
carrier = repmat(nrCarrierConfig,1,numgNBs);
for gNBIdx = 1:numgNBs
    carrier(gNBIdx).NCellID = cellIds(gNBIdx);
end
% Slot offsets of different PRS signals
prsSlotOffsets = 0:2:(2*numgNBs - 1);
prsIDs = randperm(4096,numgNBs) - 1;
% Configure PRS properties
prs = nrPRSConfig;
prs.PRSResourceSetPeriod = [10 0];
prs.PRSResourceOffset = 0;
prs.PRSResourceRepetition = 1;
prs.PRSResourceTimeGap = 1;
prs.MutingPattern1 = [];
prs.MutingPattern2 = [];
prs.NumRB = 52;
prs.RBOffset = 0;
prs.CombSize = 12;
prs.NumPRSSymbols = 12;
prs.SymbolStart = 0;
prs = repmat(prs,1,numgNBs);
for gNBIdx = 1:numgNBs
    prs(gNBIdx).PRSResourceOffset = prsSlotOffsets(gNBIdx);
    prs(gNBIdx).NPRSID = prsIDs(gNBIdx);
end

% PDSCH Configuration
pdsch = nrPDSCHConfig;
pdsch.PRBSet = 0:51;
pdsch.SymbolAllocation = [0 14];
pdsch.DMRS.NumCDMGroupsWithoutData = 1;
pdsch = repmat(pdsch,1,numgNBs);
% Path Loss Configuration
plCfg = nrPathLossConfig;
plCfg.Scenario = 'Uma';
% Generate PRS and PDSCH Resources
totSlots = nFrames*carrier(1).SlotsPerFrame;
prsGrid = cell(1,numgNBs);
dataGrid = cell(1,numgNBs);
for slotIdx = 0:totSlots-1
    [carrier(:).NSlot] = deal(slotIdx);
    [prsSym,prsInd] = deal(cell(1,numgNBs));
    for gNBIdx = 1:numgNBs
        % Create an empty resource grid spanning one slot in time domain
        slotGrid = nrResourceGrid(carrier(gNBIdx),1);

        % Generate PRS symbols and indices
        prsSym{gNBIdx} = nrPRS(carrier(gNBIdx),prs(gNBIdx));
        prsInd{gNBIdx} = nrPRSIndices(carrier(gNBIdx),prs(gNBIdx));

        % Map PRS resources to slot grid
        slotGrid(prsInd{gNBIdx}) = prsSym{gNBIdx};
        prsGrid{gNBIdx} = [prsGrid{gNBIdx} slotGrid];
    end
    % Transmit data in slots in which the PRS is not transmitted by any of
    % the gNBs (to control the hearability problem)
    for gNBIdx = 1:numgNBs
        dataSlotGrid = nrResourceGrid(carrier(gNBIdx),1);
        if all(cellfun(@isempty,prsInd))
            % Generate PDSCH indices
            [pdschInd,pdschInfo] = nrPDSCHIndices(carrier(gNBIdx),pdsch(gNBIdx));

            % Generate random data bits for transmission
            data = randi([0 1],pdschInfo.G,1);
            % Generate PDSCH symbols
            pdschSym = nrPDSCH(carrier(gNBIdx),pdsch(gNBIdx),data);

            % Generate demodulation reference signal (DM-RS) indices and symbols
            dmrsInd = nrPDSCHDMRSIndices(carrier(gNBIdx),pdsch(gNBIdx));
            dmrsSym = nrPDSCHDMRS(carrier(gNBIdx),pdsch(gNBIdx));

            % Map PDSCH and its associated DM-RS to slot grid
            dataSlotGrid(pdschInd) = pdschSym;
            dataSlotGrid(dmrsInd) = dmrsSym;
        end
        dataGrid{gNBIdx} = [dataGrid{gNBIdx} dataSlotGrid];
    end
end

% Perform OFDM modulation of PRS and data signal at each gNB
txWaveform = cell(1,numgNBs);
for waveIdx = 1:numgNBs
    carrier(waveIdx).NSlot = 0;
    txWaveform{waveIdx} = nrOFDMModulate(carrier(waveIdx),prsGrid{waveIdx} + dataGrid{waveIdx});
end

% Compute OFDM information using first carrier, assuming all carriers are
% at same sampling rate
ofdmInfo = nrOFDMInfo(carrier(1));

PRS_waveform = txWaveform{1,1};
figure %PRS自相关
grid on;
waveform_conv = xcorr(PRS_waveform);
t = linspace(-10,10,length(waveform_conv));
plot(t,abs(waveform_conv))
xlabel('时延/ms')
ylabel('归一化幅度')
% ambgfun(PRS_waveform,ofdmInfo.SampleRate,1e3,"Cut","Delay");
% ambgfun(PRS_waveform,ofdmInfo.SampleRate,1e3,"Cut","Doppler");
% Add Signal Delays and Apply Path Loss
speedOfLight = physconst('LightSpeed'); % Speed of light in m/s
sampleDelay = zeros(1,numgNBs);
radius = cell(1,numgNBs);
for gNBIdx = 1:numgNBs
   radius{gNBIdx} = rangeangle(gNBPos{gNBIdx}',UEPos');
   delay = radius{gNBIdx}/speedOfLight;                      % Delay in seconds
   sampleDelay(gNBIdx) = round(delay*ofdmInfo.SampleRate);   % Delay in samples
end
rxWaveform = zeros(length(txWaveform{1}) + max(sampleDelay),1);
rx = cell(1,numgNBs);
for gNBIdx = 1:numgNBs
    % Calculate path loss for each gNB and UE pair
    losFlag = true; % Assuming the line of sight (LOS) flag as true, as we are only considering the LOS path delays in this example
    PLdB = nrPathLoss(plCfg,fc,losFlag,gNBPos{gNBIdx}(:),UEPos(:));
    if PLdB < 0 || isnan(PLdB) || isinf(PLdB)
        error('nr5g:invalidPL',"Computed path loss (" + num2str(PLdB) + ...
            ") is invalid. Try changing the UE or gNB positions, or path loss configuration.");
    end
    PL = 10^(PLdB/10);

    % Add delay, pad, and attenuate
    rx{gNBIdx} = [zeros(sampleDelay(gNBIdx),1); txWaveform{gNBIdx}; ...
                zeros(max(sampleDelay)-sampleDelay(gNBIdx),1)]/sqrt(PL);

    % Sum waveforms from all gNBs
    rxWaveform = rxWaveform + rx{gNBIdx};
end
% TOA Estimation
cellsToBeDetected = min(3,numgNBs);
if cellsToBeDetected < 3 || cellsToBeDetected > numgNBs
    error('nr5g:InvalidNumDetgNBs',['The number of detected gNBs (' num2str(cellsToBeDetected) ...
        ') must be greater than or equal to 3 and less than or equal to the total number of gNBs (' num2str(numgNBs) ').']);
end
corr = cell(1,numgNBs);
delayEst = zeros(1,numgNBs);
maxCorr = zeros(1,numgNBs);
for gNBIdx = 1:numgNBs
    [~,mag] = nrTimingEstimate(carrier(gNBIdx),rxWaveform,prsGrid{gNBIdx});
    % Extract correlation data samples spanning about 1/14 ms for normal
    % cyclic prefix and about 1/12 ms for extended cyclic prefix (this
    % truncation is to ignore noisy side lobe peaks in correlation outcome)
    corr{gNBIdx} = mag(1:(ofdmInfo.Nfft*carrier(1).SubcarrierSpacing/15));
    % Delay estimate is at point of maximum correlation
    maxCorr(gNBIdx) = max(corr{gNBIdx});
    delayEst(gNBIdx) = find(corr{gNBIdx} == maxCorr(gNBIdx),1)-1;    
end

% Plot PRS correlation results
plotPRSCorr(corr,ofdmInfo.SampleRate);
end
    
function plotPRSCorr(corr,sr)
%   plotPRSCorr(CARRIER,CORR,SR) plots PRS-based correlation for each gNB
%   CORR, given the array of carrier-specific configuration objects CARRIER
%   and the sampling rate SR. 
    numgNBs = numel(corr);
    colors = getColors(numgNBs);
    figure;
    hold on;
    
    % Plot correlation for each gNodeB
    t = (0:length(corr{1}) - 1)/sr;
    legendstr = cell(1,2*numgNBs);
    for gNBIdx = 1:numgNBs
        plot(t,abs(corr{gNBIdx}./max(corr{gNBIdx})), ...
            'Color',colors(gNBIdx,:),'LineWidth',2);
        legendstr{gNBIdx} = sprintf('基站%d', ...
            gNBIdx);
    end

    % Plot correlation peaks
    for gNBIdx = 1:numgNBs
        c = abs(corr{gNBIdx});
        j = find(c == max(c),1);
        plot(t(j),c(j)/max(corr{gNBIdx}),'Marker','o','MarkerSize',5, ...
            'Color',colors(gNBIdx,:),'LineWidth',2);
        legendstr{numgNBs+gNBIdx} = '';
    end
    legend(legendstr);
    xlabel('时间/s');
    ylabel('归一化幅度');
    % title('PRS Correlations for All gNBs');
end

function colors = getColors(numgNBs)
%    COLORS = getColors(NUMGNBS) returns the RGB triplets COLORS for the
%    plotting purpose, given the number of gNBs NUMGNBS.

    % Form RGB triplets
    if numgNBs <= 10
        colors = [0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; ...
                  0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840; ...
                  0.9290, 0.9, 0.3; 0.9290, 0.5, 0.9; 0.6660, 0.3740, 0.2880;0, 0.4470, 0.7410];
    else
        % Generate 30 more extra RGB triplets than required. It is to skip
        % the gray and white shades as these shades are reserved for the
        % data and no transmissions in carrier grid plot in this example.
        % With this, you can get unique colors for up to 1000 gNBs.
        colors = colorcube(numgNBs+30);
        colors = colors(1:numgNBs,:);
    end
end

