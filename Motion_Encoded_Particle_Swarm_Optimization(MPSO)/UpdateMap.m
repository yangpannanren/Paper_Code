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


% 根据UAV的当前位置和测量结果，更新置信图

function [scaleFactor, newMap] = UpdateMap(currentStep,pathLength,totalMoves,dir, location, map)

    mapSize = size(map);
    mapSizeX = mapSize(1);
    mapSizeY = mapSize(2);
    
    % 目标移动方向（总路径中的totalMoves）
    move = DirToMove(dir);
    
    if totalMoves ~= 0
        moveStep = pathLength/totalMoves;
        if mod(currentStep,moveStep) == 0   % 每隔N步转移一次置信图
            tmp_map = noncircshift(map, move);
            map = tmp_map ./sum(tmp_map(:));    % 归一化
        end
    end
    
    pSensorNoDetection = ones(mapSizeX,mapSizeY); % 将没有未发现的概率初始化为1
    pSensorNoDetection(location.y,location.x) = 0; % 对于二元传感器模型，在UAV位置未被探测到的概率为零
    newMap = pSensorNoDetection .* map;   % 更新置信图
    scaleFactor = sum(newMap(:));    % 计算比例因子
    newMap = newMap ./ scaleFactor;  % 缩放更新的置信图
end