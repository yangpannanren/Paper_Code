function fig_plot_scenario(configuration,A)
% 函数功能：画定位场景构型图
% 输入参数：                            
%   configuration: 基站构型，'正方形近场','正方形远场','Y形近场','Y形远场'
%   A: 雷达和各基站的位置，第一行为雷达[x,y]，后面为基站
%   Target: 目标位置
figure
hold on;
colors = [0.4660, 0.6740, 0.1880; 0.4940, 0.1840, 0.5560; ...
          0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250;];
if(strncmp(configuration,'正',1))
    h1 = plot(A(3,1),A(3,2),'^','Color',colors(3,:),...
        'MarkerSize',11,'MarkerFace',colors(3,:),'DisplayName','基站1');
    h2 = plot(A(4,1),A(4,2),'^','Color',colors(4,:),...
        'MarkerSize',11,'MarkerFace',colors(4,:),'DisplayName','基站2');
    h3 = plot(A(2,1),A(2,2),'^','Color',colors(2,:),...
        'MarkerSize',11,'MarkerFace',colors(2,:),'DisplayName','基站3');
    h4 = plot(A(1,1),A(1,2),'s','Color','b',...
        'MarkerSize',11,'MarkerFace','b','DisplayName','雷达');
    h5 = plot(500,500,'o','Color','r',...
    'MarkerSize',11,'MarkerFace','r','DisplayName','目标');
    text(300,310,'(近场)');
    plot(5000,5000,'o','Color','r',...
    'MarkerSize',11,'MarkerFace','r','DisplayName','目标');
    text(4800,4810,'(远场)');
else
    h1 = plot(A(2,1),A(2,2),'^','Color',colors(3,:),...
        'MarkerSize',11,'MarkerFace',colors(3,:),'DisplayName','基站1');
    h2 = plot(A(4,1),A(4,2),'^','Color',colors(2,:),...
        'MarkerSize',11,'MarkerFace',colors(2,:),'DisplayName','基站3');
    h3 = plot(A(3,1),A(3,2),'^','Color',colors(1,:),...
        'MarkerSize',11,'MarkerFace',colors(1,:),'DisplayName','基站4');
    h4 = plot(A(1,1),A(1,2),'s','Color','b',...
        'MarkerSize',11,'MarkerFace','b','DisplayName','雷达');
    h5 = plot(0,500,'o','Color','r',...
    'MarkerSize',11,'MarkerFace','r','DisplayName','目标');
    text(60,510,'(近场)');
    plot(0,3000,'o','Color','r',...
    'MarkerSize',11,'MarkerFace','r','DisplayName','目标');
    text(60,3010,'(远场)');
end
xlabel('\it\fontname{Times New Roman}x/\rm\fontname{Times New Roman}m')
ylabel('\it\fontname{Times New Roman}y/\rm\fontname{Times New Roman}m')
legend([h1,h2,h3,h4,h5],'Location','northwest')
hold off;
end