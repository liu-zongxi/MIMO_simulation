%-----------------------Alamouti预编码对比-------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年4月30日09点58分-----------------%
%% 设置参数
clear;clf;
L_frame = 400; % 帧长度和帧个数
N_iter = 1500;
Nmod = 2;
M = 2^Nmod;

SNRs_dB = 0:2:20;       %信噪比
SNRs = 10.^(SNRs_dB / 10.);
N_SNR = length(SNRs_dB);

T = 2;        % Alamouti长度T       
N_M = 2;      % 列向量长度，码字长度
L_code = 64;% 码本长度
NT = 4;     % 发射和接收天线个数
NR = 1;
N_frame_bit = NT*Nmod*L_frame;  %一帧比特数
N_totoal_bit = N_frame_bit*N_iter;% 总比特数
code_book = CodebookGenerator(NT, N_M, L_code);%生成码本 

y = zeros(2,1,L_frame);
Rx = zeros(1,2,L_frame);
BERs = zeros(2, N_SNR);

fprintf('====================================================\n');
fprintf('  Precoding transmission');
fprintf('\n  %d x %d MIMO\n  %d QAM', NT,NR,M);
fprintf('\n  Simulation bits : %d', N_totoal_bit);
fprintf('\n====================================================\n');

%% 主函数
for isnr = 1:N_SNR
    SNR = SNRs(isnr);
    sigma = sqrt(NT/(2*SNR));
    n_biterror = 0;
    for iiter = 1:N_iter
        frame_origin = randi([0,1], L_frame, Nmod*N_M); % 生成随机数
        % 调制
        frame_mod = QPSKMod(frame_origin,L_frame,N_M);
        % Alamouti发送
        frame_transmit = reshape(frame_mod, N_M, 1, L_frame);
        frame_alamouti = [frame_transmit(1,1, :) -conj(frame_transmit(2,1, :));
                          frame_transmit(2,1, :) conj(frame_transmit(1,1, :))];
        % 信道和噪声
        Hiid = (randn(NR,NT)+1j*randn(NR,NT))/sqrt(2);
        arg_W = zeros(1, L_code);
        for ih = 1:L_code
            arg_W(ih) = norm(Hiid*code_book(:,:,ih),'fro');
        end
        [arg_val, arg_Index] = max(arg_W);
        HW = Hiid*code_book(:,:,arg_Index);
        for iframe = 1:L_frame
            Rx(:,:,iframe) = HW*frame_alamouti(:,:,iframe)+sigma*(randn(NR,T)+1j*randn(NR,T));
        end
        % 接收
        HW_re = [conj(HW(1)) HW(2); conj(HW(2)) -HW(1)];
        Rx_re = [Rx(1,1,:); conj(Rx(1,2,:))];
        for iframe = 1:L_frame
            y(:,:,iframe) = HW_re* Rx_re(:,:,iframe)./norm(HW)^2;
        end
        frame_re = reshape(y, L_frame, N_M);
        frame_demod = QPSKDemod(frame_re,L_frame,N_M);
        % 计算误码率
        n_biterror_tmp = sum(sum(abs(frame_demod - frame_origin)));
        n_biterror = n_biterror + n_biterror_tmp;
    end
    BERs(1, isnr) = n_biterror/(N_iter*L_frame*Nmod*N_M);
end
% 4*1 无编码
BERs(2,:) = [0.345555000000000	0.245437500000000	0.151757500000000	0.0813550000000000	0.0365075000000000	0.0133750000000000	0.00412000000000000	0.00103500000000000	0.000242500000000000	5.00000000000000e-05	2.25000000000000e-05];
semilogy(SNRs_dB,BERs(1,:),"-^", SNRs_dB,BERs(2,:),'-+');
axis([SNRs_dB([1 end]) 1e-6 1e0]);
xlabel('SNR[dB]'), ylabel('BER'); 
legend('Precoded Alamouti','No Precoded Alamouti');