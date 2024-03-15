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

% 计算与搜索路径相关的成本函数
% 输出：costP - 检测的累积概率

function costP=MyCost(position,model)
    
    if ~CheckMotion(position, model)  % 无效路径
        costP = 0;                    % 惩罚无效的路径
        return;
    else
        % 输入解
        % 注意：这里的位置实际上是一个由若干个节点组成的搜索路径
        path=PathFromMotion(position,model);
        x=path(:,1);
        y=path(:,2);

        % 输入地图
        Pmap=model.Pmap;
        N = model.n; % 路径长度
        
        % 目标移动
        targetMoves = model.targetMoves; % 目标的总动作（零表示静止）
        targetDir = model.targetDir;
        
        pNodetectionAll = 1; % 初始化完全未检测到的概率
        pDetection = zeros(N); % 初始化每个时刻的检测概率

        % 计算成本函数
        for i=1:N
            location.x = x(i) + model.xmax + 1;  % 位置被转移到[1,MAPSIZE]的范围内
            location.y = y(i) + model.ymax + 1;
            [scaleFactor,Pmap] = UpdateMap(i,N,targetMoves,targetDir,location,Pmap); % 更新概率图

            pNoDetection = scaleFactor; % 在时刻t未检测到的概率正好是缩放系数
            pDetection(i) = pNodetectionAll*(1 - pNoDetection); % 在时刻i第一次检测到的概率
            pNodetectionAll = pNodetectionAll * pNoDetection;   
        end
             
        costP = 1 - pNodetectionAll;  % 到现在为止的检测累积概率（P = 1 - R）
    end
end