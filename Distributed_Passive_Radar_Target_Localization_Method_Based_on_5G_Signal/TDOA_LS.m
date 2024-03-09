function theta = TDOA_LS(A,rho)
% 函数功能：计算TDOA_最小二乘算法定位结果
% 输入参数：                                                                 
%   A: 雷达和各基站的位置，第一行为雷达[x,y]，后面为基站
%   rho: 各基站到目标与目标到雷达的距离和
% 输出参数：
%   theta：目标定位的位置坐标[x,y]
k = 1/2*( rho.^2-vecnorm(A(2:end,:),2,2).^2 );
C = [-A(2:end,:),rho];
Y = pinv(C'*C)*C'*k;
theta=Y(1:2)';
end