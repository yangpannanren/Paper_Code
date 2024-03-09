%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TDOA定位算法LS、I2WLS、I3WLS的仿真主程序
% 地点：KC406
% 作者：UESTC-SICE-何漆龙
% 论文：《基于5G信号的分布式被动雷达目标定位方法》
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 参数设置
clear;clc;close all;
%%%%%%%%%%可调参数%%%%%%%%%%
% configuration: 基站构型，'正方形近场','正方形远场','Y形近场','Y形远场'.
configuration = 'Y形近场';
%%%%%%%%%%%%%%%%%%%%
switch configuration
    case '正方形近场'
        BS1=[0,0];BS2=[0,2];BS3=[2,0];BS4=[2,2];  %BS1为雷达坐标(km)
        Target=1e2*[5,5];
        seed = 9e4+60;
    case '正方形远场'
        BS1=[0,0];BS2=[0,2];BS3=[2,0];BS4=[2,2];
        Target=1e3*[5,5]; 
        seed = 2e5+158;
    case 'Y形近场'
        BS1=[0,0];BS2=[-2,2];BS3=[0,-2];BS4=[2,2];
        Target=1e2*[0,5];
        seed = 3e5+78;
    case 'Y形远场'
        BS1=[0,0];BS2=[-2,2];BS3=[0,-2];BS4=[2,2];
        Target=1e3*[0,3];
        seed = 1e3+461;
    otherwise
        error('configuration变量设置错误');
end
rng(seed);
std_var = logspace(-2,1.2,10); %测量噪声的标准差范围
A=1e3*[BS1;BS2;BS3;BS4]; %矩阵A包含4个初始坐标
rho = vecnorm(A(2:end,:)-Target,2,2)+vecnorm(Target,2,2); %回波路程
number=5000; %蒙特卡洛次数
RMSE=zeros(3,length(std_var));
Bias=zeros(3,length(std_var));
for j=1:length(std_var) %循环噪声标准差
    error=zeros(3,number);
    std_var1=std_var(j); %令std_var1等于当前数组的值
    for i=1:number %循环蒙特卡洛次数
        rho=rho+std_var1*randn(size(A,1)-1,1); 
        sigma=std_var1^2;
        theta1=TDOA_LS(A,rho);
        theta2=TDOA_I2WLS(A,rho,sigma);
        theta3=TDOA_I3WLS(A,rho,sigma);

        error(1,i)=norm(Target-theta1);
        error(2,i)=norm(Target-theta2);
        error(3,i)=norm(Target-theta3); 
    end
    RMSE(1,j)=(sum(error(1,:).^2)/number)^(1/2);
    RMSE(2,j)=(sum(error(2,:).^2)/number)^(1/2);
    RMSE(3,j)=(sum(error(3,:).^2)/number)^(1/2);
    
    Bias(1,j)=sum(error(1,:))/number;
    Bias(2,j)=sum(error(2,:))/number;
    Bias(3,j)=sum(error(3,:))/number;
    
    if (j==length(std_var)-1)
        % fig_plot_CDF(error) %CDF
    end
end
% RMSE(2,:)-RMSE(3,:) %直接看效果

%% 绘图
% fig_plot_scenario(configuration,A); %场景图
% fig_plot_RMSE(std_var,RMSE) %RMSE
% fig_plot_Bias(std_var,Bias) %Bias
% fig_plot_PRS_mapping; %PRS映射图
% fig_plot_PRS_match; %PRS匹配图