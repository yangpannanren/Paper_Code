function [Delay,AoA] = threeD_MUSIC(A,N_R,N,c,delat)
%% 3D-MUSIC
%========== creating seach vector=================
search_delay=(200:0.5:500)/c;
search_AoA=(0:0.005:0.5);
%========== estimation delay==============
r_delay=reshape(A(1,1,:),[],1);
R_delay=r_delay*r_delay';
[V1,~,~]=svd(R_delay);
G1=V1(:,4:N);
delay_roal=0:N-1;
for i=1:length(search_delay)
    d1=exp(-1i*2*pi*delay_roal'*delat*search_delay(i));
    p1(i)=1./abs(d1'*(G1*G1')*d1);
end
[~, peak_indices1] = findpeaks(p1, 'SortStr', 'descend');
peak1=sort(peak_indices1(1:3));
Delay=search_delay(peak1)*c;

%===============estimation AoA==================
r_AOA=A(:,1,1);
R_AOA=r_AOA*r_AOA';
[V2,~,~]=svd(R_AOA);
G1=V2(:,4:N_R);
AOA_roal=0:N_R-1;
for i=1:length(search_AoA)
    d1=exp(1i*2*pi*AOA_roal'*0.5*sin(search_AoA(i)));
    p2(i)=1./abs(d1'*(G1*G1')*d1);
end
[~, peak_indices2] = findpeaks(p2, 'SortStr', 'descend');
peak2=sort(peak_indices2(1:3));
AoA=search_AoA(peak2);
end