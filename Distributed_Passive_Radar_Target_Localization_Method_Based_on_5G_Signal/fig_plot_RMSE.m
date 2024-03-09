function fig_plot_RMSE(std_var,RMSE)
% 函数功能：画定位算法效果图-RMSE
% 输入参数：                            
%   std_var: 测量噪声的标准差范围
%   RMSE: 定位算法RMSE
figure 
semilogx(std_var,RMSE(1,:),'Color', 'b','LineWidth',2,'Marker','o','MarkerSize',8)
hold on;box on;grid on;
semilogx(std_var,RMSE(2,:),'Color', 'g','LineWidth',2,'Marker','^','MarkerSize',8);
semilogx(std_var,RMSE(3,:),'Color', 'r','LineWidth',2,'Marker','p','MarkerSize',8)
xlabel('\fontname{宋体}噪声标准差/\fontname{Times new roman}m');
ylabel('\fontname{宋体}均方根误差/\fontname{Times new roman}m');
xlim([std_var(1),std_var(end)]);
legend('LS','I2WLS','I3WLS','Location','northwest', ...
            'FontName','Times New Roman','LineWidth',1);
hold off;
end