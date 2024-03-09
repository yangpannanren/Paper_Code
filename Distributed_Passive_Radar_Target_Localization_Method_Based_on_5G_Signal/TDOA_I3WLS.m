function theta = TDOA_I3WLS(A,rho,sigma)
% 函数功能：计算TDOA_I3WLS算法定位结果
% 输入参数：
%   A: 雷达和各基站的位置，第一行为雷达[x,y]，后面为基站
%   rho: 各基站到目标与目标到雷达的距离和
%   sigma：测量噪声的方差
% 输出参数：
%   theta：目标定位的位置坐标[x,y]
m=size(A,1);
G1=[-A(2:end,:),rho];
h1=1/2*(rho.^2-vecnorm(A(2:end,:),2,2).^2); %构建矩阵h
Q=diag(ones(m-1,1)*sigma);%构建TDOA的协方差矩阵
% 初始预计值
theta0=pinv(G1'*G1)*G1'*h1;
ii=0;delta=1;
while( ii<10 && delta>1e-5 ) %初始值修正
ii=ii+1;
G3=[-2*theta0(1),-2*theta0(2),2*theta0(3)];
h3=theta0(3)^2-theta0(1)^2-theta0(2)^2;
theta3=pinv(G3'*G3)*G3'*h3;
theta0=theta0-theta3;
delta=norm(theta3(1:2)');
end
d=vecnorm(A(2:end,:)-theta0(1:2)',2,2);%根据初始估计值计算各基站到目标距离
B1=diag(d);
% 第一次wls
W1=pinv(B1*Q*B1');
theta1=pinv(G1'*W1*G1)*G1'*W1*h1;%进行第一次wls计算
cov_theta1=pinv(G1'*W1*G1);%得到theta1的协方差矩阵
%第二次wls
G2=[1,0;0,1;2*A(1,1),2*A(1,2)];%构建G'
h2=[theta1(1);theta1(2);
    theta1(1)^2+theta1(2)^2+A(1,1)^2+A(1,2)^2-theta1(3)^2];%构建h'
B2=[-1,0,0;0,-1,0;-2*theta1(1),-2*theta1(2),2*theta1(3)];%构建b'
W2=pinv(B2*cov_theta1*B2');
theta2=pinv(G2'*W2*G2)*G2'*W2*h2;
theta=theta2';%转换为（x,y）形式,得到(x0,y0)
end