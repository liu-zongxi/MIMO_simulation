%---------------------ZF 和 MMSE预编码----------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:17点18分-----------------%
%% 参数设置
clear; clf;
NT = 4;
NR = 4;     % 天线数
L_frame = 100;  %帧长度
N_iter = 1000;  % 循环次数
SNRs_dB = 3:1:30;   % 信噪比
SNRs = 10.^(SNRs_dB./10);
N_SNR = length(SNRs);
Nmod = 2;       % QPSK
N_case = 4;          % 不同类型
BERs = zeros(N_case, N_SNR);
gss = ["-kx" "-^" "-ro" "-b>" "-g<" "-m+"];   % 画图图像，注意使用双引号
y = zeros(NR,1,L_frame);
x_tilde = zeros(NR,1,L_frame);
%% 主函数
for icase = [1 2 3 4]
    gs = gss(icase);
    if icase == 1% ZF
        W_pre_formula = @(Hiid, sigma, NT) eye(size(Hiid));
        W_formula = @(Hiid, sigma, NT) inv(Hiid'*Hiid) * Hiid';
        % W_formula = @(Hiid, sigma, NT) inv(Hiid);
        beta_formula = @(Hiid, sigma, NT) 1;
    elseif icase == 2% MMSE
        W_pre_formula = @(Hiid, sigma, NT) 1;
        W_formula = @(Hiid, sigma, NT) inv(Hiid'*Hiid +sigma.^2*diag(ones(1,NT))) * Hiid';
        beta_formula = @(Hiid, sigma, NT) 1;
    elseif icase == 3
        % W_pre_formula = @(Hiid, sigma, NT) sqrt(NT/trace(inv(Hiid)*inv(Hiid)'))*inv(Hiid);
        W_pre_formula = @(Hiid, sigma, NT) sqrt(NT/trace((Hiid' * inv(Hiid*Hiid'))*(Hiid' * inv(Hiid*Hiid'))'))*Hiid' * inv(Hiid*Hiid');
        W_formula = @(Hiid, sigma, NT) eye(size(Hiid));
        beta_formula = @(Hiid, sigma, NT) sqrt(NT/trace((Hiid' * inv(Hiid*Hiid'))*(Hiid' * inv(Hiid*Hiid'))'));
    elseif icase == 4
        W_pre_formula = @(Hiid, sigma, NT) sqrt(NT/trace((Hiid'*inv(Hiid*Hiid'+sigma.^2*eye(NR,NR)))*(Hiid'*inv(Hiid*Hiid'+sigma.^2*eye(NR,NR)))'))*Hiid'*inv(Hiid*Hiid'+sigma.^2*eye(NR,NR));
        W_formula = @(Hiid, sigma, NT) eye(size(Hiid));
        beta_formula = @(Hiid, sigma, NT) sqrt(NT/trace((Hiid'*inv(Hiid*Hiid'+sigma.^2*eye(NR,NR))*(Hiid'*inv(Hiid*Hiid'+sigma.^2*eye(NR,NR)))')));
    end
    for isnr = 1:N_SNR
        SNR = SNRs(isnr);
        n_biterror = 0;
        for iiter = 1:N_iter
            % 生成数据
            frame_origin = randi([0,1],L_frame,Nmod*NT);
            % QPSK调制
            frame_mod=QPSKMod(frame_origin,L_frame, NT);
            % 调整为标准形式
            frame_transmit = reshape(frame_mod, NT, 1, L_frame);
            % 生成信道，SIMO有NR个信道
            Hiid = (randn(NR,NT)+1j*randn(NR,NT))./sqrt(2);
            %  AWGN噪声
            sigma = sqrt(NT/(2*SNR));
            noise = sigma*(randn(NR, 1) + 1j*randn(NR, 1));
            % 预编码矩阵
            W_pre = W_pre_formula(Hiid, sigma, NT);
            % 接收信号
            beta = beta_formula(Hiid, sigma, NT);
            for iframe = 1:L_frame
%                 if icase == 1
%                     y(:,:,iframe) = (Hiid*frame_transmit(:,:,iframe)+noise)/beta;
%                 end
                y(:,:,iframe) = (Hiid * W_pre * frame_transmit(:,:,iframe)+noise)/beta;
            end
            % 信号检测
            W = W_formula(Hiid, sigma, NT);
            for iframe = 1:L_frame
                x_tilde(:,:,iframe) = W * y(:,:,iframe);
            end
            % 解调
            frame_re = reshape(x_tilde, L_frame, NT);
            frame_demod = QPSKDemod(frame_re,L_frame,NT);
            % 计算误码率
            n_biterror_tmp = sum(sum(abs(frame_demod - frame_origin)));
            n_biterror = n_biterror + n_biterror_tmp;
        end
        BERs(icase, isnr) = n_biterror/(N_iter*L_frame*Nmod*2);
    end
end
%% 画图
semilogy(SNRs_dB,BERs(1,:),gss(1), SNRs_dB,BERs(2,:),gss(2), SNRs_dB,BERs(3,:),gss(3), SNRs_dB,BERs(4,:),gss(4));
axis([SNRs_dB([1 end]) 1e-6 1e0])
xlabel('SNR[dB]'), ylabel('BER'); 
legend('ZF','MMSE', "Pre-ZF", "Pre-MMSE");