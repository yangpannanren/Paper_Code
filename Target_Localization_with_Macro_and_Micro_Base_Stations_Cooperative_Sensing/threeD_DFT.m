function [estiamte_target1,estiamte_target2,estiamte_target3] = threeD_DFT(A,c,delta_f,Nr,N,MBS_location,MiBS_location)
A=squeeze(A(:,1,:));
fft_result_dim1 = fft(A,[],1);
DD = ifft(fft_result_dim1,[],2);
% x=1:N;
% y=1:Nr;
% [X,Y]=meshgrid(x,y);
% mesh(X,Y,abs(DD))
% figure(2)
% plot(sum(abs(DD),2))
% findpeaks(sum(abs(DD),2))
% figure(3)
% plot(sum(abs(DD),1))
% findpeaks(sum(abs(DD),1))
%=================== AoA searching==============
[peaks1, peak_indices1] = findpeaks(sum(abs(DD),2));
[~, sorted_peak_indices1] = sort(peaks1, 'descend');
dim_AoA = sort(peak_indices1(sorted_peak_indices1(1:3)));
%=================delay searching================
[peaks2, peak_indices2] = findpeaks(sum(abs(DD),1));
[~, sorted_peak_indices2] = sort(peaks2, 'descend');
dim_delay =sort( peak_indices2(sorted_peak_indices2(1:3)));

%============ transform the index to estiamtion============
dimA=zeros(1,3);
dimC=zeros(1,3);
for i=1:3
    dimA(i)=asin(2*(dim_AoA(i)-1)/Nr);
end
for i=1:3
    dimC(i)=(dim_delay(i)-1)*c/(delta_f*N);
end

%============ AoA localization====================
distance_MBS_MiBS=sqrt((MBS_location(1)-MiBS_location(1))^2+(MBS_location(2)-MiBS_location(2))^2);
% target1
Tar_1=(dimC(1)^2-distance_MBS_MiBS^2)/(2*(dimC(1)-distance_MBS_MiBS*cos(dimA(1))));
estiamte_target1=[MBS_location(1)+Tar_1*cos(dimA(1)),MBS_location(2)+Tar_1*sin(dimA(1))];
% target
Tar_2=(dimC(2)^2-distance_MBS_MiBS^2)/(2*(dimC(2)-distance_MBS_MiBS*cos(dimA(2))));
estiamte_target2=[MBS_location(1)+Tar_2*cos(dimA(2)),MBS_location(2)+Tar_2*sin(dimA(2))];
% target
Tar_3=(dimC(3)^2-distance_MBS_MiBS^2)/(2*(dimC(3)-distance_MBS_MiBS*cos(dimA(3))));
estiamte_target3=[MBS_location(1)+Tar_3*cos(dimA(3)),MBS_location(2)+Tar_3*sin(dimA(3))];
end