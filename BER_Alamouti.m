%------------------Alamouti分集编码的BER性能----------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:22点28分-----------------%
%% 设置参数
clear;clf;
SNRs_dB = 0:2:20;   % 信噪比
SNRs = 10.^(SNRs_dB./10);
N_SNR = length(SNRs_dB);    % 信噪比个数
N_iter = 1000;      % 迭代次数
Nmod = 2;           % 调制阶数
L_frame = 100;      % 帧长度
Niter = 1000;       % 迭代次数
N_case = 6;          % 不同类型
BERs = zeros(N_case, N_SNR);
gss = ["-kx" "-^" "-ro" "-b>" "-g<" "-m+"];   % 画图图像，注意使用双引号
%% 主函数
for icase = 1:N_case
    if icase == 1   % AWGN 信道
        NR = 1;
        NT = 1;
        power_T = 2;
        y_formula = @(Hiid, frame, noise) frame + noise;
    elseif icase == 2   % SISO瑞利衰落信道,可以用MRC代替
        NR = 1;
        NT = 1;
        power_T = 2;
        % y_formula = @(Hiid, frame, noise) frame + noise./Hiid;
        y_formula = @(Hiid, frame, noise) sum(conj(Hiid).*(Hiid.*frame + noise),2);
    elseif icase == 3   % 1*2 MRC瑞利衰落信道
        NR = 2;
        NT = 1;
        power_T = 2;
        % y_formula = @(Hiid, frame, noise) frame + noise./Hiid;
        y_formula = @(Hiid, frame, noise) sum(conj(Hiid).*(Hiid.*frame + noise), 2);
    elseif icase == 4   % 2*1 Alamouti编码
        NR = 1;
        NT = 2;
        power_T = 1;
        y_formula = @(Hiid, frame, noise) frame + [conj(Hiid(:,1)).* noise(:,1)+Hiid(:,2).*conj(noise(:,2)) conj(Hiid(:,2)).*noise(:,2)-Hiid(:, 1).*conj(noise(:,2))]./(sum(abs(Hiid.^2), 2));
    elseif icase == 5
        NR = 2;
        NT = 2;
        power_T = 1;
        y_formula = @(Hiid, frame, noise) frame + [conj(Hiid(:,1)).* noise(:,1)+conj(Hiid(:,2)).*noise(:,2)+Hiid(:,3).*conj(noise(:,3))+Hiid(:, 4).*conj(noise(:,4)) conj(Hiid(:,3)).*noise(:, 1)+conj(Hiid(:,4)).*noise(:, 2)-Hiid(:,1).*conj(noise(:,3))-Hiid(:,2).*conj(noise(:,4))]./(sum(abs(Hiid.^2), 2));
    elseif icase == 6
        NR = 1;
        NT = 4;
        power_T = 1;
        y_formula = @(Hiid, frame, noise) frame + [conj(Hiid(:,1)).* noise(:,1)+conj(Hiid(:,2)).*noise(:,2)+conj(Hiid(:,3)).*noise(:,3)+conj(Hiid(:,4)).*noise(:,4) conj(Hiid(:,2)).*noise(:,1)-conj(Hiid(:, 1)).*noise(:,2)-conj(Hiid(:, 4)).*noise(:,3)+conj(Hiid(:, 3)).*noise(:,4) conj(Hiid(:,3)).*noise(:,1)+conj(Hiid(:, 4)).*noise(:,2)-conj(Hiid(:, 1)).*noise(:,3)-conj(Hiid(:, 2)).*noise(:,4) conj(Hiid(:,4)).*noise(:,1)-conj(Hiid(:, 3)).*noise(:,2)+conj(Hiid(:, 2)).*noise(:,3)-conj(Hiid(:, 1)).*noise(:,4)]./(sum(abs(Hiid.^2), 2));
    end
    gs = gss(icase);
    for isnr = 1:N_SNR
        SNR = SNRs(isnr);
        n_biterror = 0;
        for iiter = 1:N_iter
            % 生成数据,不管几根天线，都发一组数据，因为这是分集而不是MIMO
            frame_origin = randi([0,1],L_frame,Nmod*NT);
            % QPSK调制
            frame_mod=QPSKMod(frame_origin,L_frame, NT);
            % 生成信道，SIMO有NR个信道
            Hiid = (randn(L_frame,NR*NT)+1j*randn(L_frame,NR*NT))/sqrt(2);
            %  AWGN噪声
            sigma = sqrt(1/(2*power_T*SNR));
            noise = sigma*(randn(L_frame, NR*NT) + 1j*randn(L_frame, NR*NT));
            y = y_formula(Hiid, frame_mod, noise);
            % 解调
            frame_demod = QPSKDemod(y,L_frame,NT);
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
legend('AWGN信道','SISO瑞利衰落信道', '1发2收MCR方案', '2发1收Alamouti方案', '2发2收Alamouti方案')
xlabel('信噪比Eb/N0')
ylabel('误比特率（BER）')
title('2发2收Alamouti方案在瑞利衰落信道下的性能')