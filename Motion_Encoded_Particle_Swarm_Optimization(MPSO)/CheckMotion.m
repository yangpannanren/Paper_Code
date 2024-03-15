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

% 检查编码后的运动是否形成了有效路径
%
function valid = CheckMotion(position, model)
    
    % 加载模型参数
    n=model.n;
    xs = model.xs; % 初始搜索位置
    ys = model.ys;
    path = zeros(n,2);  % 路径由n个节点组成，每个节点为(x,y)
    currentNode = [xs ys]; % 当前节点
    valid = true; % 有效
    
    % 从动作到路径的解码
    for i=1:n
        motion = position(i,:);
        nextMove = MotionDecode(motion); % 将编码的运动从8个数字序列中转换为动作
        nextNode = currentNode + nextMove;
        % 在地图边界之外
        if nextNode(1) > model.xmax || nextNode(1) < model.xmin...
            || nextNode(2) > model.ymax || nextNode(2) < model.ymin
            valid = false;
            return
        end            
        path(i,:) = nextNode;
        currentNode = nextNode;
    end
   
    % 检查重复的行
    [u,I,J] = unique(path, 'rows', 'first'); % 将 path 中的每一行视为单个实体，并按排序顺序返回 path 中的唯一行
    if size(u,1) < size(path,1)
        valid = false;
        return
    end

end