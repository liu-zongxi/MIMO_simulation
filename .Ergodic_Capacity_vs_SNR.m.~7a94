%------------------MIMO信道容量和信噪比的影响---------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:21点43分-----------------%
% Ergodic_Capacity_vs_SNR.m
%% 参数设置
clear;
SNRs_dB=[0:5:20];       % 信噪比
SNRs=10.^(SNRs_dB/10.); 
N_SNR = length(SNRs);
N_iter=1000;            % 迭代次数
NTRs = [1 1; 1 2; 2 1; 4 2; 4 4];   % 天线个数，每一行是一种情况
Ncase = size(NTRs, 1);              % 测试不同情况的个数
C = zeros(Ncase, N_SNR);
%% 主函数
for iicase = 1:Ncase
    % 初始化天线个数
    NT = NTRs(iicase, 1);
    NR = NTRs(iicase, 2);
    rank = min(NT, NR);
    I = eye(rank);
    for iiiter = 1:N_iter
        H = sqrt(0.5)*(randn(NR,NT)+1j*randn(NR,NT));
        % 为了和单位矩阵保持一致
        if NR>=NT
            HH = H'*H; 
        else
            HH = H*H'; 
        end
        for iiSNR = 1:N_SNR
            C(iicase,iiSNR) = C(iicase,iiSNR)+log2(real(det(I+SNRs(iiSNR)/NT*HH)));
        end
    end
end
C = C/N_iter;
plot(SNRs_dB,C(1,:),'b-o', SNRs_dB,C(2,:),'b-<', SNRs_dB,C(3,:),'b-s',SNRs_dB,C(4,:),'b->', SNRs_dB,C(5,:),'b-^');
xlabel('SNR[dB]'); 
ylabel('bps/Hz'); 
set(gca,'fontsize',10); 
grid on
s1='{\it N_T}=1,{\it N_R}=1'; 
s2='{\it N_T}=1,{\it N_R}=2'; 
s3='{\it N_T}=2,{\it N_R}=1'; 
s4='{\it N_T}=2,{\it N_R}=2'; 
s5='{\it N_T}=4,{\it N_R}=4';
legend(s1,s2,s3,s4,s5)
title('未知CSI时的MIMO信道遍历容量')