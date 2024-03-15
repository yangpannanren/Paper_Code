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

% 绘制搜索地图

function PlotModel(model)
    
    MAP_SIZE = model.MAPSIZE;
    x = 1:1:MAP_SIZE; y = 1:1:MAP_SIZE; % x & y从0到20
    [X,Y] = meshgrid(x,y); % 复制x & y以创建一个矩形网格(X,Y)
    clf('reset');
    h = pcolor(X,Y,model.Pmap);
    hold on;

    % 起始位置
    plot(model.xs+MAP_SIZE/2+0.5,model.ys+MAP_SIZE/2+0.5,'w-o',...
        'MarkerSize',3,...
        'MarkerFaceColor','w',...
        'LineWidth',1);
    xlabel('x (cell)');
    ylabel('y (cell)');
    
    % 设置属性
    set(h,'linewidth',0.1);
    cb = colorbar; 
    cb.Ruler.Exponent = -3;
    gcaP=get(gca,'position');
    cbP=get(cb,'Position');
    cbP(3)=cbP(3)/2; % 改变色条宽度
    set(cb,'Position',cbP)
    set(gca,'position',gcaP)
    set(gcf,'position',[300,100,350,250]); % 设置地图尺寸
end