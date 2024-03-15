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

%
% 创建随机路径（解）
% 

function position=CreateRandomSolution(model)

    n=model.n;  % 加载路径节点的数量
    startNode = [model.xs model.ys];
    path = zeros(n,2);  % 路径初始化
    MaxIt = 100; % 重置路径前的最大试验迭代次数

    motions = [ 1       0;
                0.7071  0.7071;
                0       1;
               -0.7071  0.7071;
               -1       0;
               -0.7071 -0.7071;
                0      -1;
                0.7071 -0.7071];
             
    should_restart = true;
    
    % 重复进行，直到生成一个有效的路径
    while should_restart
        should_restart = false;
        for i = 1:n
            path(i,:) = startNode;
        end
        position = zeros(n,2); % 运动初始化
        currentNode = startNode;
        for i=1:n
            motion = motions(randi([1 length(motions)]),:);
            invalidFlag = true;
            it = 0;
            while (invalidFlag && it < MaxIt)
                nextMove = MotionDecode(motion);
                nextNode = currentNode + nextMove;
                invalidFlag = false;

                % 将路径限制在地图内
                % x方向超出 -> 移回
                if nextNode(1) > model.xmax
                    motion = motions(5,:);
                    invalidFlag = true;
                    it = it + 1;
                elseif nextNode(1) < model.xmin
                    motion = motions(1,:);
                    invalidFlag = true;
                    it = it + 1;
                    
                % y方向超出
                elseif nextNode(2) > model.ymax
                    motion = motions(7,:);
                    invalidFlag = true;
                    it = it + 1;
                elseif nextNode(2) < model.ymin
                    motion = motions(3,:);
                    invalidFlag = true;
                    it = it + 1;
                else
                    % 检查路径中的重复节点
                    for j=1:length(path)
                        if isequal(nextNode,path(j,:))
                            motion = motions(randi([1 length(motions)]),:);
                            invalidFlag = true;
                            it = it + 1;
                            break;
                        end
                    end 
                end
            end
               
            % 重新开始整个路径
            if (it >= MaxIt)
                should_restart = true;
                break;
            else    % 路径 ok
                path(i,:) = nextNode;
                currentNode = nextNode;
                position(i,:) = motion;
            end
        end
    end
end