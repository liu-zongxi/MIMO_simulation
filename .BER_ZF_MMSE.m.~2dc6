%------------------ZF和MMSE检测算法----------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:17点18分-----------------%
%% 参数设置
NT = 4;
NR = 4;     % 天线数
L_frame = 100;  %帧长度
N_iter = 1000;  % 循环次数
SNRs_dB = 3:1:20;   % 信噪比
SNRs = 10.^(SNRs_dB./10);
N_SNR = length(SNRs);
Nmod = 2;       % QPSK
N_case = 2;          % 不同类型
BERs = zeros(N_case, N_SNR);
gss = ["-kx" "-^" "-ro" "-b>" "-g<" "-m+"];   % 画图图像，注意使用双引号
%% 主函数
for icase = 1:N_case
    gs = gss(icase);
    if icase == 1
        W_formula = @(Hiid, sigma, NT) Hiid'*inv(Hiid*Hiid');
    elseif icase == 2
        W_formula = @(Hiid, sigma, NT) Hiid'*inv(Hiid*Hiid'+2*sigma.^2*diag(ones(1,NT)));
    end
    for isnr = 1:N_SNR
        SNR = SNRs(isnr);
        n_biterror = 0;
        for iiter = 1:N_iter
            % 生成数据
            frame_origin = randi([0,1],L_frame,Nmod*NT);
            % QPSK调制
            frame_mod=QPSKMod(frame_origin,L_frame, NT);
            % 生成信道，SIMO有NR个信道
            Hiid = (randn(NR,NT)+1j*randn(NR,NT))./sqrt(2);
            %  AWGN噪声
            sigma = sqrt(NT/(2*SNR));
            noise = sigma*(randn(L_frame, NR) + 1j*randn(L_frame, NR));
            % 接收信号
            y = frame_mod*Hiid+noise;
            % 信号检测
            W = W_formula(Hiid, sigma, NT);
            x_tilde = y*W;
            % 解调
            frame_demod = QPSKDemod(x_tilde,L_frame,NT);
            % 计算误码率
            n_biterror_tmp = sum(sum(abs(frame_demod - frame_origin)))
            n_biterror = n_biterror + n_biterror_tmp;
        end
        BERs(icase, isnr) = n_biterror/(N_iter*L_frame*Nmod*2);
        semilogy(SNRs_dB,BERs(icase,:),gs);
        hold on;
        axis([SNRs_dB([1 end]) 1e-6 1e0])
    end
end