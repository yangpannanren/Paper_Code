% 添加路径
addpath('Routines')
addpath('Algorithms')

% -----------------------------------------------------------
% 说明
% -----------------------------------------------------------

% 一个实现。
% 源、噪声和多路径信道的位置都是随机的。

%% ----------------------------------------------------------
%   场景参数
% -----------------------------------------------------------
% 所有基站的天线数量
Nantennas = 70*ones(4,1);

% 载波频率 [Hz]
carrierFreq = 7e9;

% 带宽 [Hz]
bandwidth = 30e6;

% 地图的宽度
mapWidth = 100;

% 基站的位置
baseStations = [45,45; 45,-45; -45,45; -45,-45];

% 校准误差
phaseStdDev = 0; % 相位标准差
gainStdDevDB = 0; % 增益标准差[dB]

% 多径参数
clusterDecayTime = 34e-9; %簇衰减时间
rayDecayTime = 29e-9; %射线衰减时间
clusterInterTime = 17e-9; %簇间时间
rayInterTime = 5e-9; %射线间时间
angleStdDev = 26/180*pi; %角度标准差

%% ----------------------------------------------------------
%   算法参数
% -----------------------------------------------------------

% 定位技术选择（只选一个）
%   - 11 仪器变量估算器（IV）
%   - 21 基于平方距离的最小二乘估计（SR-LS）
%   - 31 Stansfield 估计器
%   - 41 直接位置确定（DPD）
%   - 42 DiSouL

algorithm = 42;

% 网格技术的起始分辨率 [m]
startCellSize = 5;

% 每个天线的LOS脉冲能量与噪声频谱密度之比 [dB]
EdivN0 = 10;

% 过采样因子
oversamplingFactor = 3;

% 生成新信号、新噪声和多径信道？
randomize = true;


%% ----------------------------------------------------------
%    开始模拟
% -----------------------------------------------------------
switch algorithm
    case 11
        temp = 'Instrumental Variables estimator (IV)';
    case 21
        temp = 'Squared-Range-based Least Squares';
    case 31
        temp = 'Stansfield estimator';
    case 41
        temp = 'DPD';
    case 42
        temp = 'DiSouL';
end
disp('% -----------------------------------------------------------')
disp(['%   ', temp, ' results'])
disp('% -----------------------------------------------------------')

% 生成随机信号和传感器位置
Nstations = length(Nantennas);
if randomize
    % 生成均匀环形阵列位置和
    antLocs = cell(Nstations,1);
    fraunhoferDist = zeros(Nstations,1);
    for i=1:Nstations
        angleSep = 2*pi/Nantennas(i); %天线之间的角度间隔
        radiusArray = 3e8/carrierFreq/4/sin(pi/Nantennas(i)); %阵列半径
        orienation = rand*2*pi; %随机的方位角
        % 天线的位置
        antLocs{i} = radiusArray*[cos(angleSep*(1:Nantennas(i)).'+orienation),sin(angleSep*(1:Nantennas(i)).'+orienation)];
        antLocs{i} = ones(Nantennas(i),1)*baseStations(i,:) +antLocs{i};
        % 弗洛恩费尔距离
        fraunhoferDist(i) = 2*(2*radiusArray)^2/3e8*carrierFreq; 
    end
    
    % 发射源的位置
    mobileUser = .5*mapWidth*(rand(1,2)-.5);
    % countTrials 是一个用于控制循环次数的计数器。
    % 它被用于确保生成的移动用户位置 mobileUser 不在基站的Fraunhofer距离范围内。
    % 小于于Fraunhofer距离称为天线的近场，大于为远场。天线特性通常在远场测量。
    % 如果无法在规定次数内找到适当的移动用户位置，就会抛出一个错误。
    % 这样可以确保生成的移动用户位置符合模拟特定条件的要求。
    countTrials = 0; 
    while any(sqrt(sum((baseStations-ones(Nstations,1)*mobileUser).^2,2))<fraunhoferDist)
        countTrials=countTrials+1;
        mobileUser = .5*mapWidth*(rand(1,2)-.5);
        if countTrials>100
            error('The arrays are too large!!')
        end
    end

    % 生成多径信道
    [trueTOAs, trueAOAs, trueAmplis ] = generateChannels( mobileUser, antLocs, ones(Nstations,1), clusterDecayTime, rayDecayTime, clusterInterTime, rayInterTime, angleStdDev );

    % 计算噪声功率
    noisePower = 1/10^(EdivN0/10); %信号功率归一化为1

    % 对所有天线应用MF
    % 第一项表示信号传播的时间，考虑了信号从发射到接收所需的时间以及过采样因子等因素。
    % 第二项表示簇衰减时间、带宽和过采样因子的乘积，考虑了多径信道的影响，需要观察足够长的时间来捕获多径效应。
    % 第三项表示过采样因子的影响。
    % 将这些项相加，然后向上取整，得到了观察到足够长时间以确保接收信号能量变得可以忽略的样本数。
    Nsamples = ceil(sqrt(2)*mapWidth/3e8*bandwidth*oversamplingFactor +20*clusterDecayTime*bandwidth*oversamplingFactor +2*oversamplingFactor);
    % % % % % 不是很懂
    
    [receivedSignals, receivedSignalsAfterMF, startIndex ] = generateUncalibratedSignals( antLocs, carrierFreq, bandwidth, oversamplingFactor, noisePower, Nsamples, trueTOAs, trueAOAs, trueAmplis, phaseStdDev, gainStdDevDB );

end

% 通过阈值法MF估计TOA（正偏）
[TOAbyThresMF,~] = globalThresMF(receivedSignalsAfterMF,bandwidth,oversamplingFactor,startIndex,noisePower);

% 通过找到首次超过阈值的时间估计TOA
[TOAbyThres1stCross,snapshots] = globalFirstCrossThresMF(receivedSignalsAfterMF,bandwidth,oversamplingFactor,startIndex,noisePower);

% 丢弃没有检测到信号的BS数据
select = ~cellfun(@isempty,snapshots);
snapshots = snapshots(select);
activeTOAbyThresMF = TOAbyThresMF(select);
activeTOAbyThres1stCross = TOAbyThres1stCross(select);
activeBSLocs = baseStations(select,:);
activeantLocs = antLocs(select);

% 估计AOA
% cellfun(@(x)x(1)/pi*180,trueAOAs)
AOAs = beamforming( activeantLocs, snapshots, carrierFreq);

% 运行定位算法
switch algorithm
    case 11
        mobileUserEstim = IV( activeBSLocs, AOAs );
    case 21
        mobileUserEstim = SR_LS( activeBSLocs, activeTOAbyThresMF );
    case 31
        mobileUserEstim = Stansfield( activeBSLocs, activeTOAbyThresMF, AOAs );
    case 41 % DPD
        mobileUserEstim = DPD( antLocs, receivedSignals, carrierFreq, bandwidth, oversamplingFactor, mapWidth, startCellSize, 5 );
    case 42 % AOA-based Direct Positioning Assisted by TOA
        [mobileUserEstim,losBS] = DiSouL( activeantLocs, snapshots, noisePower, startCellSize, carrierFreq, bandwidth, activeTOAbyThresMF, mapWidth );
end


%% ----------------------------------------------------------
%   输出结果
% -----------------------------------------------------------

disp(['The mobile user is located at [', num2str(mobileUser,'%1.1f '),'] m.'])
if isempty(mobileUserEstim)
    disp('The mobile user''s location could not be estimated.')
else
    rMSE = norm(mobileUserEstim-mobileUser);
    disp(['The estimated mobile user''s location is [', num2str(mobileUserEstim,'%1.1f '),'] m.'])
    disp(['The localization error is ', num2str(rMSE,3),' m.'])
end

% 生成2D图形
opengl software
clf('reset');
set(gcf,'Color','white');
hold on
% 添加以估计TOA为半径的圆
for l=1:length(activeTOAbyThresMF)
    if isnumeric(activeTOAbyThresMF(l))
        [x,y] = circle(activeBSLocs(l,:),3e8*activeTOAbyThresMF(l),100);
        patch(x,y,[.8,0,0],'EdgeColor','none','FaceAlpha',1/sum(~isnan(activeTOAbyThresMF)));
    end
end
% 添加表示近场的圆
for l=1:Nstations
    [x,y] = circle(baseStations(l,:),fraunhoferDist(l),100);
    patch(x,y,[.8,0,0],'EdgeColor','none','FaceColor','white');
end
% 添加表示数据快照中到达方向的真实方向线
snapshotsAmplis = findSnapshotAmplis( bandwidth, trueTOAs(select), trueAmplis(select), activeTOAbyThres1stCross );
for l=1:Nstations
    if isnumeric(TOAbyThres1stCross(l))
        select = abs(snapshotsAmplis{l}).^2*Nantennas(l) > noisePower;
        selectedAOAs = trueAOAs{l}(select);
        for i=1:length(selectedAOAs)
            plot(baseStations(l,1)+[0,2*mapWidth*cos(selectedAOAs(i))],baseStations(l,2)+[0,2*mapWidth*sin(selectedAOAs(i))],'-b')
        end
    end
end
hBS = scatter(baseStations(:,1),baseStations(:,2),100,'fill','o','MarkerEdgeColor','none','MarkerFaceColor','g');
hMU = scatter(mobileUser(1),mobileUser(2),150,'fill','s','MarkerEdgeColor','k','MarkerFaceColor','w');
if ~isempty(mobileUserEstim)
    hMUestim = scatter(mobileUserEstim(1),mobileUserEstim(2),70,'fill','s','MarkerEdgeColor','k','MarkerFaceColor','c');
    legend([hBS,hMU,hMUestim],{'Nodes','Source','Source estimate'})
else
    legend([hBS,hMU],{'Nodes','Source'})
end
set(gca,'box','on');
set(gca,'layer','top');
hold off
xlim([-mapWidth/2,mapWidth/2]);
ylim([-mapWidth/2,mapWidth/2]);
xlabel('x-axis [m]');
ylabel('y-axis [m]');

% print('../results/Fig1','-dpng')