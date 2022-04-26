%------------------瑞利信道下MRC分集的性能----------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:22点25分-----------------%
%% 参数设置
clear;clf;
L_frame= 130;   % 帧长度
N_packet = 4000;% 包个数
Nmod = 2;       % 调制阶数
M = 2^Nmod;     % 
SNRs_dB = [0:2:20]; % 信噪比
SNRs = 10.^(SNRs_dB/10);
N_SNR = length(SNRs);
NTRs = [1 1; 1 2; 1 4];   % 天线个数，每一行是一种情况
N_case = size(NTRs, 1);              % 测试不同情况的个数
BERs = zeros(N_case, N_SNR);
BERs_y = zeros(N_case, N_SNR);
gss = ["-kx" "-^" "-ro"];   % 画图图像，注意使用双引号
%% 主函数
for icase = 1:N_case
    NT = NTRs(icase, 1);
    NR = NTRs(icase, 2);
    gs = gss(icase);
    for isnr = 1:N_SNR
        n_biterror = 0;
        SNR = SNRs(isnr);
        sigma = sqrt(0.5/SNR);  % 噪信比，实数虚数各占一半能量，用于噪声幅度
        for ipacket = 1:N_packet
            % 生成数据
            frame_origin = randi([0,1],L_frame,NT*Nmod);
            % QPSK调制
            frame_mod=QPSKMod(frame_origin,L_frame, NT);
            % 生成信道，SIMO有NR个信道
            % 何为瑞利信道？就是乘性的一个瑞利衰落信道
            Hiid = (randn(L_frame,NR)+1j*randn(L_frame,NR))/sqrt(2);
            noise = sigma*(randn(L_frame,NR) +1j * randn(L_frame,NR));
            % 接收信号
            y = Hiid .* frame_mod + noise;  % 这就是公式y=SNR*x+z,但归一化为系数在z之前
            % 进行MRC
            W_mrc = conj(Hiid);
            y_mrc = sum(W_mrc.*y,2);
            % 解调
            frame_demod = QPSKDemod(y_mrc,L_frame,NT);
            % 计算误码率
            n_biterror_tmp = sum(sum(abs(frame_demod - frame_origin)));
            n_biterror = n_biterror + n_biterror_tmp;
        end
        BERs(icase, isnr) = n_biterror/(N_packet*L_frame*Nmod);
        semilogy(SNRs_dB,BERs(icase,:),gs);
        hold on;
        axis([SNRs_dB([1 end]) 1e-6 1e0])
    end
end
title('BER perfoemancde of MRC Scheme');
xlabel('SNR[dB]');
ylabel('BER') 
grid on;
set(gca,'fontsize',9)
legend('SISO','MRC (Tx:1,Rx:2)','MRC (Tx:1,Rx:4)')
