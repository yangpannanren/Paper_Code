function fig_plot_CDF(error)
% 函数功能：画定位算法效果图-CDF
% 输入参数：                            
%   error: 定位算法误差
figure
hold on;box on;grid on;
h1 = cdfplot(error(1,:));
h2 = cdfplot(error(2,:));
h3 = cdfplot(error(3,:));
set(h1,'LineStyle', '--', 'Color', 'b','LineWidth',2)
set(h2,'LineStyle', '-.', 'Color', 'g','LineWidth',2)
set(h3,'LineStyle', '-', 'Color', 'r','LineWidth',2)
xlabel('\fontname{宋体}定位误差/\fontname{Times new roman}m');
ylabel('累积分布函数','FontName','宋体');
delete(get(gca,'title'));
legend('LS','I2WLS','I3WLS','Location','southeast', ...
    'FontName','Times New Roman','LineWidth',1);
hold off
end