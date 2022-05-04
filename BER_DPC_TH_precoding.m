%------------------------DPC和TH编码-------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年5月4日10点38分-----------------%
%% 参数设置
clear;clf;
L_frame =100;
N_packet = 2000;
Nmod = 2;
NT = 4;
N_user = 10;
N_user_ative = 4;
I=eye(NT,NT);

SNRs_dB = 0:2:20;
SNRs = 10.^(SNRs_dB/10);
N_SNR = length(SNRs);

BERs = zeros(2, N_SNR);

Ncase = 2;
gss = ["-kx" "-b^" "-ro" "-m+"];   % 画图图像，注意使用双引号
%% 主函数
for icase = 1:Ncase
    gs = gss(icase);
    for isnr = 1:N_SNR
        SNR = SNRs(isnr);
        n_biterror = 0;
        sigma = sqrt(NT/(2*SNR));
        for ipacket = 1:N_packet
            % 生成数据
            frame_origin = randi([0,1],L_frame,NT*Nmod);
            % QPSK调制
            frame_mod=QPSKMod(frame_origin,L_frame, NT);
            frame_reshape = reshape(frame_mod,NT,L_frame);
            % 信道
            H = (randn(N_user,NT)+1j*randn(N_user,NT))/sqrt(2);
            % 选出四个最好的信道
            Types = nchoosek([1:N_user],N_user_ative);
            for itype = 1:size(Types, 1)
                H_sel = H(Types(itype,:)',:);
                [Q_H, R_H] = qr(H_sel);
                L_min(itype) = min(diag(R_H));
            end
            % 最大的最小值就是最好的四个
            [val_max_L_min,index_max_L_min] = max(L_min);
            H_selected = H(Types(index_max_L_min,:)',:);
            % Matlab没有LQ分解，要这样曲线救国
            [Q_temp,R_temp] = qr(H_selected');
            L=R_temp';
            Q=Q_temp';
            % 预编码
            frame_prepre = frame_reshape;
            if icase == 1
                for ipre = 2:4
                    frame_prepre(ipre,:) = frame_prepre(ipre,:) - L(ipre,1:ipre-1)/L(ipre,ipre)*frame_prepre(1:ipre-1,:);
                end
            else
                for ipre = 2:4
                    frame_prepre(ipre,:) = ModTH(frame_prepre(ipre,:) - L(ipre,1:ipre-1)/L(ipre,ipre)*frame_prepre(1:ipre-1,:), 2);
                end
            end
            frame_pre = Q'*frame_prepre;
            % 信道
            noise = sigma*(randn(N_user_ative,L_frame)+1j*randn(N_user_ative,L_frame));
            y = H_selected*frame_pre+noise;
            % 接收
            frame_re = inv(diag(diag(L)))*y;
            frame_recieved = reshape(frame_re,NT*L_frame,1);
            if icase==2 % in the case of TH precoding
                frame_recieved = ModTH(frame_recieved,2);
            end
            frame_pre_demod = reshape(frame_recieved,L_frame,NT);
            frame_demod = QPSKDemod(frame_pre_demod,L_frame,NT);
            % 计算误码率
            n_biterror_tmp = sum(sum(abs(frame_demod - frame_origin)))
            n_biterror = n_biterror + n_biterror_tmp;
        end
        BERs(icase, isnr) = n_biterror/(N_packet*L_frame*Nmod);
    end
    semilogy(SNRs_dB,BERs(icase,:),gs);
    hold on;
end