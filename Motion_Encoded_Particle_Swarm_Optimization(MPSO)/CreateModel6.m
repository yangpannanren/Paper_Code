%_________________________________________________________________________%
%  Motion-encoded Partical Swarm Optimization (MPSO) source codes demo 1.0%
%                                                                         %
%  Developed in MATLAB 2020a                                              %
%                                                                         %
%  Author and programmer: Manh Duong Phung                                %
%                                                                         %
%         e-Mail: duongpm@gmail.com                                       %
%                 duongpm@vnu.edu.vn                                      %
%                                                                         %
%       Homepage: https://uet.vnu.edu.vn/~duongpm/                        %
%                                                                         %
%   Main paper: Manh Duong Phung, Quang Phuc Ha                           %
%               "Motion-encoded particle swarm optimization for moving    %
%               target search using UAVs",                                %
%               Applied soft computing , in press,                        %
%               DOI: https://doi.org/10.1016/j.asoc.2020.106705           %
%                                                                         %
%_________________________________________________________________________%

% 创建具有初始置信度的搜索地图
% 场景6与场景5类似，但起始位置在潜在区域的下方，目标向东北方向移动。
function model=CreateModel6()

    %% 创建一个网格地图
    MAP_SIZE = 40; % 地图大小
    x = 1:1:MAP_SIZE;
    y = 1:1:MAP_SIZE; 
    [X,Y] = meshgrid(x,y); % 复制x和y以创建一个矩形网格(X,Y)

    % 生成概率图
    mu1 = [15 8];  % 第一个分布的均值（潜在的目标位置）
    Sigma1 = MAP_SIZE*[0.1 0;0 0.1]; % 分布的协方差
    F1 = mvnpdf([X(:) Y(:)],mu1,Sigma1); % F1中各点的多元正态分布的概率密度函数(pdf)
    F1 = reshape(F1,length(y),length(x)); % 将F1转换为一个矩阵
    F1 = F1/sum(F1(:)); % 归一化
    
    mu2 = [8 15];  % 第二个分布/信息的均值（潜在的目标位置）
    Sigma2 = MAP_SIZE*[0.08 0;0 0.08]; % 分布的协方差
    F2 = mvnpdf([X(:) Y(:)],mu2,Sigma2); % F2中各点的多元正态分布的概率密度函数(pdf)
    F2 = reshape(F2,length(y),length(x)); % 将F2转换为一个矩阵  
    F2 = F2/sum(F2(:)); % 归一化
    
    Pmap = 0.5*(F1 + F2); % 对有关目标的两个信息源的概率进行标记
    
    %% 设置
    % 地图限制
    xmin= -floor(MAP_SIZE/2);
    xmax= floor(MAP_SIZE/2);
    
    ymin= -floor(MAP_SIZE/2);
    ymax= floor(MAP_SIZE/2);
    
    % 初始搜索位置
    xs=-12;
    ys=-12;
    
    % 路径节点的数量（不包括起始位置（起始节点））
    n=20;
    
    % 动作范围
    MRANGE = 4;
    
    %% 将地图和搜索参数纳入一个模型中
    model.xs=xs;
    model.ys=ys;
    model.Pmap=Pmap;
    model.n=n;
    model.xmin=xmin;
    model.xmax=xmax;
    model.ymin=ymin;
    model.ymax=ymax;
    model.MRANGE = MRANGE;
    model.MAPSIZE = MAP_SIZE;
    model.X = X;
    model.Y = Y;
    model.targetMoves = 10; % 必须能被路径长度整除（例如，mod(N,move)=0）
    model.targetDir = 'NE';  % 目标向东北移动

end