%------------------已知和未知CSI的信道容量对比---------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:22点07分-----------------%
%OL_CL_Comparison.m
%% 设置参数
clear;clc;
SNRs_dB=[0:5:20];       % 信噪比
SNRs=10.^(SNRs_dB/10.);
N_SNR = length(SNRs);   % 信噪比长度
N_iter=1000;            % 迭代次数
NT=4;                   % 4*4矩阵
NR=4; 
rank=min(NT,NR);        % MIMO的秩
I = eye(rank);          % 单位矩阵
rho=0.2;                % 相关系数
Rtx=[1      rho     rho^2   rho^3;  % 发射相关矩阵
     rho     1      rho    rho^2;
     rho^2   rho     1       rho;  
     rho^3   rho^2   rho     1];
Rrx=[1      rho     rho^2   rho^3;  % 接受相关矩阵
    rho     1       rho     rho^2;
    rho^2   rho     1       rho; 
    rho^3   rho^2   rho     1];
C_OL=zeros(1,length(SNRs_dB));
C_CL=zeros(1,length(SNRs_dB));
%% 主函数    
for iiter=1:N_iter
   Hiid = sqrt(1/2)*(randn(NT,NR) + 1j*randn(NT,NR));   % 生成一个独立同分布H
   H = Rrx^(1/2)*Hiid*Rtx^(1/2);  % 窄带信道
   sigma = svd(H'*H);
   for i=1:N_SNR
      %random channel generation
      C_OL(i) = C_OL(i) + log2(det(I+SNRs(i)*(H'*H)/NT));
      % Gamma = Water_Pouring(sigma,SNRs(i),NT);
      Gamma = WaterFilling(H, rank, SNRs(i), NT);
      C_CL(i) = C_CL(i)+log2(det(I+SNRs(i)/NT*diag(Gamma)*diag(sigma)));
   end
end
C_OL = real(C_OL)/N_iter;  
C_CL = real(C_CL)/N_iter;
figure, plot(SNRs_dB, C_OL,'-o', SNRs_dB, C_CL,'-<');
xlabel('SNR [dB]');
ylabel('bps/Hz'); 
set(gca,'fontsize',10);
legend('Channel Unknown','Channel Known'); 
title('开环和闭环MIMO信道容量')
grid on