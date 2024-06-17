clear all
%% general parameters
Nc=256;  % number of subcarrier
M_sym=1; % number of symbols
delta_f1=30e3; % subcarriers spacing of MBS
delta_f2=120e3; % subcarriers spacing of MiBS
c_0=3e8; % light speed
I=3; % number of  target
N_transmit=32; % number of transmit antenna
N_receive=32; % number of receive antenna
a=0.5; % the antenna space/ wavelength
MC=10; % monte carlo times
Q=delta_f2/delta_f1; % the ratio of subcarrier spacing
%% target and BS parameter
target_1=[200,30]; % the target 1 location
target_2=[250,60]; % the target 2 location
target_3=[300,80]; % the target 3 location
MBS_location=[0,0]; % MBS location
MiBS_location=[300,0]; % MiBS location

%% Delay of targets
%============ Total delay of target 1 =========
Range1=sqrt((target_1(1)-MiBS_location(1))^2+(target_1(2)-MiBS_location(2))^2)+sqrt((target_1(1)-MBS_location(1))^2+(target_1(2)-MBS_location(2))^2);
delay1=Range1/c_0; % delay of target 1 to MiBS and MBS
%============ Total delay of target 2 =========
Range2=sqrt((target_2(1)-MiBS_location(1))^2+(target_2(2)-MiBS_location(2))^2)+sqrt((target_2(1)-MBS_location(1))^2+(target_2(2)-MBS_location(2))^2);
delay2=Range2/c_0; % delay of target 2 to MiBS and MBS
%============ Total delay of target 3 =========
Range3=sqrt((target_3(1)-MiBS_location(1))^2+(target_3(2)-MiBS_location(2))^2)+sqrt((target_3(1)-MBS_location(1))^2+(target_3(2)-MBS_location(2))^2);
delay3=Range3/c_0; % delay of target 3 to MiBS and MBS

%% AoA and AoD of targets
%=========== target 1=================
AoA_1=atan(target_1(2)/target_1(1)); % target to MBS
AoD_1=atan(target_1(2)/(MiBS_location(1)-target_1(1))); % target to MiBS
%=========== target 2=================
AoA_2=atan(target_2(2)/target_2(1)); % target to MBS
AoD_2=atan(target_2(2)/(MiBS_location(1)-target_2(1))); % target to MiBS
%=========== target 3=================
AoA_3=atan(target_3(2)/target_3(1)); % target to MBS
AoD_3=atan(target_3(2)/(MiBS_location(1)-target_3(1))); % target to MiBS

%% 3D matrix from MBS
%============ 3D matrix of target1 from MBS ============
three_dimentional_target1_MBS=zeros(N_receive,N_transmit,Nc); % creating a 3D matrix
AoA_target1=zeros(N_receive,1);
AoD_target1=zeros(1,N_transmit);
for k=1:N_receive
    AoA_target1(k,1)=exp(1i*2*pi*a*(k-1)*sin(AoA_1));
end
for p=1:N_transmit
    AoD_target1(1,p)=exp(1i*2*pi*a*(p-1)*sin(AoD_1));
end
AoA_AoD_2D_target1=AoA_target1*AoD_target1; % creating a AoA-AoD matrix
for n=1:Nc
    three_dimentional_target1_MBS(:,:,n)=AoA_AoD_2D_target1.*exp(-1i*2*pi*(n-1)*delta_f2*delay1);
end

%============ 3D matrix of target2 from MBS ============
three_dimentional_target2_MBS=zeros(N_receive,N_transmit,Nc); % creating a 3D matrix
AoA_target2=zeros(N_receive,1);
AoD_target2=zeros(1,N_transmit);
for k=1:N_receive
    AoA_target2(k,1)=exp(1i*2*pi*a*(k-1)*sin(AoA_2));
end
for p=1:N_transmit
    AoD_target2(1,p)=exp(1i*2*pi*a*(p-1)*sin(AoD_2));
end
AoA_AoD_2D_target2=AoA_target2*AoD_target2; % creating a AoA-AoD matrix
for n=1:Nc
    three_dimentional_target2_MBS(:,:,n)=AoA_AoD_2D_target2.*exp(-1i*2*pi*(n-1)*delta_f2*delay2);
end

%============ 3D matrix of target3 from MBS ============
three_dimentional_target3_MBS=zeros(N_receive,N_transmit,Nc); % creating a 3D matrix
AoA_target3=zeros(N_receive,1);
AoD_target3=zeros(1,N_transmit);
for k=1:N_receive
    AoA_target3(k,1)=exp(1i*2*pi*a*(k-1)*sin(AoA_3));
end
for p=1:N_transmit
    AoD_target3(1,p)=exp(1i*2*pi*a*(p-1)*sin(AoD_3));
end
AoA_AoD_2D_target3=AoA_target3*AoD_target3; % creating a AoA-AoD matrix
for n=1:Nc
    three_dimentional_target3_MBS(:,:,n)=AoA_AoD_2D_target3.*exp(-1i*2*pi*(n-1)*delta_f2*delay3);
end

%% 3D matrix from MiBS
%============ 3D matrix of target1 from MiBS ============
three_dimentional_target1_MiBS=zeros(N_receive,N_transmit,Nc); % creating a 3D matrix
AoA_target1=zeros(N_receive,1);
AoD_target1=zeros(1,N_transmit);
for k=1:N_receive
    AoA_target1(k,1)=exp(1i*2*pi*a*(k-1)*sin(AoD_1));
end
for p=1:N_transmit
    AoD_target1(1,p)=exp(1i*2*pi*a*(p-1)*sin(AoA_1));
end
AoA_AoD_2D_target1=AoA_target1*AoD_target1; % creating a AoA-AoD matrix
for n=1:Nc
    three_dimentional_target1_MiBS(:,:,n)=AoA_AoD_2D_target1.*exp(-1i*2*pi*(n-1)*delta_f1*delay1);
end

%============ 3D matrix of target2 from MiBS ============
three_dimentional_target2_MiBS=zeros(N_receive,N_transmit,Nc); % creating a 3D matrix
AoA_target2=zeros(N_receive,1);
AoD_target2=zeros(1,N_transmit);
for k=1:N_receive
    AoA_target2(k,1)=exp(1i*2*pi*a*(k-1)*sin(AoD_2));
end
for p=1:N_transmit
    AoD_target2(1,p)=exp(1i*2*pi*a*(p-1)*sin(AoA_2));
end
AoA_AoD_2D_target2=AoA_target2*AoD_target2; % creating a AoA-AoD matrix
for n=1:Nc
    three_dimentional_target2_MiBS(:,:,n)=AoA_AoD_2D_target2.*exp(-1i*2*pi*(n-1)*delta_f1*delay2);
end

%============ 3D matrix of target3 from MiBS ============
three_dimentional_target3_MiBS=zeros(N_receive,N_transmit,Nc); % creating a 3D matrix
AoA_target3=zeros(N_receive,1);
AoD_target3=zeros(1,N_transmit);
for k=1:N_receive
    AoA_target3(k,1)=exp(1i*2*pi*a*(k-1)*sin(AoD_3));
end
for p=1:N_transmit
    AoD_target3(1,p)=exp(1i*2*pi*a*(p-1)*sin(AoA_3));
end
AoA_AoD_2D_target3=AoA_target3*AoD_target3; % creating a AoA-AoD matrix
for n=1:Nc
    three_dimentional_target3_MiBS(:,:,n)=AoA_AoD_2D_target3.*exp(-1i*2*pi*(n-1)*delta_f1*delay3);
end

%% RMSE
RMSE_target1=[];
for snr=-20:5:20
    RMSE1=0;
    tic
    for m=1:MC
        %% Add AWGN noise 
        %================== MBS side=====================
        three_dimentional_all_MBS=(three_dimentional_target1_MBS+three_dimentional_target2_MBS+three_dimentional_target3_MBS);
        three_dimentional_all_MBS_noise=awgn(three_dimentional_all_MBS,snr);
        %================== MiBS side====================
        three_dimentional_all_MiBS=(three_dimentional_target1_MiBS+three_dimentional_target2_MiBS+three_dimentional_target3_MiBS);
        three_dimentional_all_MiBS_noise=awgn(three_dimentional_all_MiBS,snr+30);
        %% fusion
        %========= step1 :  transposing the 3D-MiBS matrix========
        New_MiBS_noise=permute(three_dimentional_all_MiBS_noise, [2, 1, 3]);
        
        %========= fusion =================
        D=zeros(N_receive,N_transmit,Q*Nc);
        D(:,:,1:Nc)=D(:,:,1:Nc)+New_MiBS_noise(:,:,:); % add the processed MiBS matrix
        for n=1:Nc
               D(:,:,Q*n-(Q-1))=(D(:,:,Q*n-(Q-1))+three_dimentional_all_MBS_noise(:,:,n));
        end
%         %% sensing
%         %===========fusion sensing===============
%         [Total_delay,Total_AoA]=threeD_MUSIC(D,N_receive,Nc,c_0,delta_f1);
% 
%         %% calculating the location
%         %=========== fusion ==================
%         [tarF_1,tarF_2,tarF_3]=AoA_location(MBS_location,MiBS_location,Total_delay,Total_AoA);
        [tarF_1,tarF_2,tarF_3]=grid_3D_DFT(D,N_receive,Q*Nc,delta_f1,c_0,MBS_location,MiBS_location);
%          [tarF_1,tarF_2,tarF_3]=threeD_DFT(D,c_0,delta_f1,N_receive,Q*Nc,MBS_location,MiBS_location);
          RMSE1=RMSE1+((tarF_3(1)-target_3(1))^2+(tarF_3(2)-target_3(2))^2+ ...
            (tarF_2(1)-target_2(1))^2+(tarF_2(2)-target_2(2))^2+ ...
            (tarF_1(1)-target_1(1))^2+(tarF_1(2)-target_1(2))^2)/3;
    end
    toc
    RMSE_target1=[RMSE_target1,sqrt(RMSE1/MC)];
end
