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


% 从编码的运动中创建一个搜索路径
%
function path=PathFromMotion(position,model)

    n=model.n;
    xs = model.xs;
    ys = model.ys;
    path = zeros(n,2);  % 路径初始化
    currentNode = [xs ys];
    for i=1:n
        motion = position(i,:);
        nextMove = MotionDecode(motion);
        nextNode = currentNode + nextMove;
        
        % 将路径限制在地图内
        % x 方向
        if nextNode(1) > model.xmax
            nextNode(1) = model.xmax;
        elseif nextNode(1) < model.xmin
            nextNode(1) = model.xmin;
        end
        % y 方向
        if nextNode(2) > model.ymax
            nextNode(2) = model.ymax;
        elseif nextNode(2) < model.ymin
            nextNode(2) = model.ymin;    
        end
        path(i,:) = currentNode;
        currentNode = nextNode;
    end
end