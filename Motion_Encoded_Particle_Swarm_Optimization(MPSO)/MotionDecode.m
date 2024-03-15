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

% 将编码的运动从8个数字序列中转换为动作
% 8个数字编码一个矢量的大小，其角度为：
% NW, N, NE, W, E, SW, S, SE
% 1   2   3  4  5   6  7  8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function move = MotionDecode(motion)
    
    % 找出和向量的角度
    angle = atan2(motion(2), motion(1));
    
    % 将角度映射到其相应的八度空间（0，45，90，135...）
    octant = mod(round( 8 * angle / (2*pi) + 8 ),8) + 1;
    moveArray = [1  0;
                 1  1;
                 0  1;
                -1  1;
                -1  0;
                -1 -1;
                 0 -1;
                 1 -1];
     % 将角度映射到一个移动上
     move = moveArray(octant,:);
end