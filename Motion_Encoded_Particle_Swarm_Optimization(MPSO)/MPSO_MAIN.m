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


% 主程序：运动编码粒子群优化（MPSO）
%
% 寻找一条找到目标的概率最大的路径
% 

clc;
clear;
close all;

%% 创建搜索场景

% 创建搜索地图和参数
model = CreateModel1(); % 场景1
% model = CreateModel2(); % 场景2
% model = CreateModel3(); % 场景3
% model = CreateModel4(); % 场景4
% model = CreateModel5(); % 场景5
% model = CreateModel6(); % 场景6

CostFunction=@(x) MyCost(x,model);  % 成本函数

nVar = model.n;       % 决策变量的数量=PSO的搜索维度=运动的数量

VarSize=[nVar 2];   % 决策变量矩阵的大小

VarMin=-model.MRANGE;           % 粒子的下限（变量）
VarMax = model.MRANGE;          % 粒子的上界 

%% PSO 参数

MaxIt=100;          % 最大迭代次数

nPop=1000;          % 种群大小 (粒子群大小)

w=1;                % 惯性权重
wdamp=0.98;         % 惯性权重阻尼比
c1=2.5;             % 个体学习因子
c2=2.5;             % 全局学习因子

alpha= 2;
VelMax=alpha*(VarMax-VarMin);    % 粒子最大速度
VelMin=-VelMax;                  % 粒子最小速度

%% 初始化

% 创建一个空的粒子结构
empty_particle.Position=[];
empty_particle.Velocity=[];
empty_particle.Cost=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];

% 初始化全局最优
GlobalBest.Cost = -1; % 最大化问题

% 创建一个空的粒子矩阵，每个粒子是一个解（搜索路径）
particle=repmat(empty_particle,nPop,1);

% 初始化循环
for i=1:nPop
    
    % 初始化粒子位置
    particle(i).Position=CreateRandomSolution(model);
    
    % 初始化粒子速度
    particle(i).Velocity=zeros(VarSize);
    
    % 评估
    costP = CostFunction(particle(i).Position);
    particle(i).Cost= costP;
    
    % 更新个体最优
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    
    % 更新全局最优
    if particle(i).Best.Cost>GlobalBest.Cost
        GlobalBest=particle(i).Best;
    end
    
end

% 在每个迭代中保持最优成本函数值的数组
BestCost=zeros(MaxIt,1);

%% PSO主循环
for it=1:MaxIt
    for i=1:nPop
                
        % 更新粒子速度
        particle(i).Velocity = w*particle(i).Velocity ...
            + c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            + c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % 更新粒子速度范围
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % 更新粒子位置
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % 更新粒子位置范围
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        % 评估
        costP = CostFunction(particle(i).Position);
        particle(i).Cost = costP;

        % 更新个人最优
        if particle(i).Cost > particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            % 更新全局最优
            if particle(i).Best.Cost > GlobalBest.Cost
                GlobalBest=particle(i).Best;
            end
            
        end
      
    end
    
    % 更新找到的最优成本函数
    BestCost(it)=GlobalBest.Cost;
    
    % 惯性权重阻尼比
    w=w*wdamp;

    % 显示迭代信息
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
end

%% 结果

% 根据目标移动情况更新地图
targetMoves = model.targetMoves; % 目标移动的数量（零表示静态）
moveDir = DirToMove(model.targetDir); % 目标运动的方向
moveArr = targetMoves*moveDir;
updatedMap = noncircshift(model.Pmap, moveArr);
newModel = model;
newModel.Pmap = updatedMap;

% 绘制解
figure(1);
path=PathFromMotion(GlobalBest.Position,model); % 从运动空间转换到笛卡尔空间
PlotSolution(path,newModel);  

% 绘制迭代过程中的最优成本函数图
figure(2);
plot(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
