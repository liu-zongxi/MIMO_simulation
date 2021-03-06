%------------------MU-MIMO的信道反转方式---------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年5月3日15点50分-----------------%
%% 参数设置
clear;clf;
L_frame = 100;      % 一个帧里多少符号
N_packet = 2000;    % 一个包里多少帧
Nmod = 2;           % 调制阶数
NT = 4;             % 发射天线数
N_user = 20;        % 用户数
N_user_active = NT;  % 当前使用的用户数
I = eye(NT,NT);     
Ncase = 4;

SNRs_dB = 0:2:20;
SNRs = 10.^(SNRs_dB/10);
N_SNR = length(SNRs);
BERs = zeros(Ncase, N_SNR);
gss = ["-kx" "-b^" "-ro" "-m+"];   % 画图图像，注意使用双引号
%% 主函数
for icase = 1:Ncase
    gs = gss(icase);
    if icase ==1
        W_formula = @(H, sigma, I) H'*inv(H*H');
        flag_select = 0;
    elseif icase ==2
        W_formula = @(H, sigma, I) H'*inv(H*H'+sigma*I);
        flag_select = 0;
    elseif icase ==3
        W_formula = @(H, sigma, I) H'*inv(H*H');
        flag_select = 1;
    elseif icase ==4
        W_formula = @(H, sigma, I) H'*inv(H*H'+sigma*I);
        flag_select = 1;
    end 
    for isnr = 1:N_SNR
        SNR = SNRs(isnr);
        n_biterror = 0;
        sigma = sqrt(NT/(2*SNR));
        for ipacket = 1:N_packet
            % 生成信号
            frame_origin = randi([0,1],L_frame,N_user_active*Nmod);
            frame_mod = QPSKMod(frame_origin,L_frame, N_user_active);
            frame_reshape = reshape(frame_mod, N_user_active, L_frame);
            % 生成信道
            if flag_select == 0
                H_DL = (randn(N_user_active,NT)+1j*randn(N_user_active,NT))/sqrt(2);
            elseif flag_select ==1
                H = zeros(N_user, NT);
                Channel_norm = zeros(1, N_user);
                for iuser = 1:N_user
                    H(iuser,:) = (randn(1,NT)+1j*randn(1,NT))/sqrt(2);
                    Channel_norm(iuser)=norm(H(iuser,:));
                end
                [Ch_norm,Index]=sort(Channel_norm,'descend');
                H_DL = H(Index(1:N_user_active),:);
            end
            W = W_formula(H_DL, sigma, I);
            beta = sqrt(NT/trace(W*W'));
            W_pre = beta*W;
            frame_trnasmit = W_pre*frame_reshape;
            noise = sigma*(randn(N_user_active,L_frame)+1j*randn(N_user_active,L_frame));
            y = H_DL * frame_trnasmit + noise;
            % 接收
            frame_re = y/beta;
            frame_pre_demod = reshape(frame_re, L_frame, N_user_active);
            frame_demod = QPSKDemod(frame_pre_demod,L_frame,N_user_active);
            % 计算误码率
            n_biterror_tmp = sum(sum(abs(frame_demod - frame_origin)))
            n_biterror = n_biterror + n_biterror_tmp;
        end
        BERs(icase, isnr) = n_biterror/(N_packet*L_frame*Nmod);
    end
    semilogy(SNRs_dB,BERs(icase,:),gs);
    hold on;
end
%% 画图
title('MU-MIMO BER inv');
xlabel('SNR[dB]');
ylabel('BER') 
grid on;
set(gca,'fontsize',9)
legend('ZF','MMSE','ZF with select', 'MMSE with select')
