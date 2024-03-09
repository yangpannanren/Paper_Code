function fig_plot_PRS_mapping()
% 函数功能：画5G_PRS映射图
% 参考链接：https://ww2.mathworks.cn/help/5g/ug/5g-new-radio-prs.html
% Set carrier parameters
carrier = nrCarrierConfig;
% Set parameters related to PRS slot configuration
prs = nrPRSConfig;
prs.PRSResourceSetPeriod = [8 0];
prs.PRSResourceOffset = [0 4];
prs.PRSResourceRepetition = 2;
prs.PRSResourceTimeGap = 1;
% numSlots = 43;                                       % Consider 43 slots to compare the plot with Figure 1
% plotTitle = 'PRS Slot Configuration';
% plotGrid(carrier,prs,numSlots,'SlotFill',plotTitle); % Slot numbers on the x-axis of the MATLAB plot are 1-based
prs.MutingPattern1 = [1 1];                          % Use [] to disable the muting bit pattern option-1
prs.MutingBitRepetition = 2;
prs.MutingPattern2 = [1 0];                          % Use [] to disable the muting bit pattern option-2
carrier.NSizeGrid = 1;
prs = nrPRSConfig;
prs.NumPRSSymbols = 12;
prs.SymbolStart = 1;
prs.NumRB = 1;
prs.RBOffset = 0;
prs.CombSize = 4;
prs.REOffset = 2;
numSlots = 1;                                      % Consider one slot to highlight the RE pattern in the MATLAB plot
plotTitle = 'PRS映射结构';
plotGrid(carrier,prs,numSlots,'REFill',plotTitle); % OFDM symbol and subcarrier numbers in the MATLAB plot are 1-basedgrid on;
end

function plotGrid(carrier,prs,numSlots,flag,imTitle)
    [transmittedSlots,mutedSlots] = getTransmittedAndMutedSlots(carrier,prs,numSlots); % Of size numSlots-by-numRes
    numRes = size(transmittedSlots,2);
    rng(37);                  % Set RNG state for repeatability
    tempColors = reshape(rand((numRes-2)*3,1),[],3);
    map = [1 1 1;...          % White color
           0 0 0.6667;...     % Blue color
           1 0 0;...          % Red color
           tempColors];
    
    resScalings = 1:numRes+1; % 1 for white color, 2 for blue color, and 3 for red color, and so on
    prsREGrid = zeros(carrier.NSizeGrid*12,carrier.SymbolsPerSlot*numSlots);
    prsSlotGrid = resScalings(1)*ones(carrier.NSizeGrid*12,numSlots);                  % For white background

    figure();
    % Plot PRS slot grid
    image(abs(prsSlotGrid));

    % Apply colormap to image
    colormap(map);

    for slotIdx = 0:numSlots-1
        if strcmpi(flag,'REFill')
            carrier.NSlot = slotIdx;
            ind = nrPRSIndices(carrier,prs,'OutputResourceFormat','cell');
            sym = nrPRS(carrier,prs,'OutputResourceFormat','cell');
            slotGrid = nrResourceGrid(carrier);
            for resIdx = 1:numel(ind)
                slotGrid(ind{resIdx}) = resScalings(resIdx+1)*abs(sym{resIdx});
            end
            prsREGrid(:,(1:carrier.SymbolsPerSlot)+carrier.SymbolsPerSlot*slotIdx) = slotGrid;

            % Replace all zero values of carrier grid with proper scaling
            % for white background
            prsREGrid(prsREGrid == 0) = resScalings(1);

            % Plot carrier grid
            image(abs(prsREGrid));
            axis xy;

            % Apply colormap to image
            colormap(map);

            % Add labels to x-axis and y-axis
            xlabel('OFDM符号');
            ylabel('子载波');

            % Generate lines
            L = line(ones(numRes),ones(numRes),'LineWidth',8);
            % Index color map and associate selected colors with lines
            set(L,{'color'},mat2cell(map(resScalings(2:end),:),ones(1,numRes),3));
            % Create legend
            legendNames = cell(1,numRes);
            for resIdx = 1:numRes
                legendNames{resIdx} = ['PRS resource #' num2str(resIdx)];
            end
            % legend(legendNames{:});
        else % 'SlotFill'
            axis xy;
            hold on;
            for resIdx = 1:numRes
                ismuted = mutedSlots(slotIdx+1,resIdx);
                isTransmitted = transmittedSlots(slotIdx+1,resIdx);
                temp = [slotIdx+1 slotIdx+1];
                color = map(resIdx+1,:);
                if isTransmitted
                    hT(resIdx) = patch([temp(1)-0.5 temp(1)-0.5 temp(end)+0.5 temp(end)+0.5],[1 624 624 1],...
                        color,'LineStyle','none'); %#ok<AGROW>
                end
                if ismuted
                    hM(resIdx) = patch([temp(1)-0.5 temp(1)-0.5 temp(end)+0.5 temp(end)+0.5],[1 624 624 1],...
                        color,'LineStyle','none','FaceAlpha',0.5); %#ok<AGROW>
                end
            end
        end
    end

    % Add title to image
    % title(imTitle);
    if strcmpi(flag,'SlotFill')
        % Add labels to x-axis and y-axis
        xlabel('Slots');
        ylabel('Subcarriers');

        % Create legend
        legendNames = cell(1,numRes);
        legendNamesTxRes = cell(1,numRes);
        legendNamesMutedRes = cell(1,numRes);
        for resIdx = 1:numRes
            % legendNames{resIdx} = ['PRS resource #' num2str(resIdx)];
            legendNamesTxRes{resIdx} = ['Transmitted instance of PRS resource #' num2str(resIdx)];
            legendNamesMutedRes{resIdx} = ['Muted instance of PRS resource #' num2str(resIdx)];
        end
        if sum(mutedSlots(:)) > 0
            legend([hT hM],[legendNamesTxRes legendNamesMutedRes]);
        else
            legend(hT,legendNames);
        end
    end
end

function [transmittedSlots,mutedSlots] = getTransmittedAndMutedSlots(carrier,prs,numSlots)
%   [TRANSMITTEDSLOTS,MUTEDSLOTS] = getTransmittedAndMutedSlots(CARRIER,PRS,NUMSLOTS)
%   returns logical arrays to give information about transmitted slots
%   TRANSMITTEDSLOTS and muted slots MUTEDSLOTS for all PRS resources by
%   considering these inputs:
%
%   CARRIER  - Carrier specific configuration object
%   PRS      - Positioning reference signal configuration object
%   NUMSLOTS - Number of slots for which the output is returned

    % Take copy of input PRS configuration to retain input muting
    % configuration unchanged for further processing
    prs1 = prs;

    % Extract PRS configuration without considering muting aspects
    prs1.MutingPattern1 = [];
    prs1.MutingBitRepetition = 2;
    prs1.MutingPattern2 = [];
    numResOffVal = numel(prs1.PRSResourceOffset);
    numSymStartVal = numel(prs1.SymbolStart);
    numPRSSymVal = numel(prs1.NumPRSSymbols);
    numREOffVal = numel(prs1.REOffset);
    numNPRSIDVal = numel(prs1.NPRSID);

    % Calculate the number of PRS resources configured in a PRS
    % resource set
    numRes = max([numResOffVal, numSymStartVal, numPRSSymVal, ...
        numREOffVal, numNPRSIDVal]);
    PRSPresenceWithOutMuting = zeros(numSlots,numRes);
    for slotIdx = 0:numSlots-1
        carrier.NSlot = slotIdx;
        [~,PRSPresence] = nr5g.internal.validateAndSchedulePRS(carrier,prs1);
        PRSPresenceWithOutMuting(slotIdx+1,:) = PRSPresence;
    end

    % Consider PRS configuration with muting aspects
    prs1.MutingPattern1 = prs.MutingPattern1;
    prs1.MutingBitRepetition = prs.MutingBitRepetition;
    prs1.MutingPattern2 = prs.MutingPattern2;
    PRSPresenceWithMuting = zeros(numSlots,numRes);
    for slotIdx = 0:numSlots-1
        carrier.NSlot = slotIdx;
        [~,PRSPresence] = nr5g.internal.validateAndSchedulePRS(carrier,prs1);
        PRSPresenceWithMuting(slotIdx+1,:) = PRSPresence;
    end

    % Identify transmitted and muted slots based on PRS scheduling
    transmittedSlots = PRSPresenceWithMuting;                                 % Of size numSlots-by-numRes
    mutedSlots = PRSPresenceWithOutMuting ~= PRSPresenceWithMuting;           % Of size numSlots-by-numRes
end