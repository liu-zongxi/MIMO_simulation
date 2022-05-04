%------------------MU-MIMO的对角块化---------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年5月3日15点50分-----------------%
%% 设置
clear;
clf;
L_frame = 100;
N_packet = 20000;
Nmod = 2;
NT = 4;
NR = 2;
N_user = 2;

gss = ["-kx" "-^" "-ro" "-b>" "-g<" "-m+"];
SNRs_dB = 0:2:20;
SNRs = 10.^(SNRs_dB/10);
N_SNR = length(SNRs_dB);
BERs = zeros(5, N_SNR);
%% 主函数
for isnr = 1:N_SNR
    SNR = SNRs(isnr);
    n_biterror = 0;
    sigma = sqrt(NT/(2*SNR));
    for ipacket = 1:N_packet
        % 生成信号
        frame_origin = randi([0,1],L_frame,NT*Nmod);
        frame_mod = QPSKMod(frame_origin,L_frame, NT);
        frame_reshape = reshape(frame_mod, NT, L_frame);
        % 信道
        noise = sigma* (randn(NR,L_frame)+1j*randn(NR,L_frame));
        H1 = (randn(NR,NT)+1j*randn(NR,NT))/sqrt(2);
        H2 = (randn(NR,NT)+1j*randn(NR,NT))/sqrt(2);
        H = (randn(NR,NT, N_user)+1j*randn(NR,NT, N_user))/sqrt(2);
        % 对角化处理
        for i = 1:N_user

        end
        [U1,S1,V1] = svd(H1);
        W2 = V1(:,3:4);
        [U2,S2,V2] = svd(H2); 
        W1 = V2(:,3:4);
        frame_transmit = W1*frame_reshape(1:2,:) + W2*frame_reshape(3:4,:);
        y1 = H1*frame_transmit+noise;
        y2 = H2*frame_transmit+noise;
        HV1 = H1*W1;
        EQ1 = HV1'*inv(HV1*HV1'); % Equalizer for the 1st user
        HV2 = H2*W2;
        EQ2 = HV2'*inv(HV2*HV2'); % Equalizer for the 2nd user
        y = [EQ1*y1; EQ2*y2];
        frame_pre_demod = reshape(y,L_frame,NT);
        frame_demod = QPSKDemod(frame_pre_demod,L_frame,NT);
        % 计算误码率
        n_biterror_tmp = sum(sum(abs(frame_demod - frame_origin)))
        n_biterror = n_biterror + n_biterror_tmp;
    end
    BERs(1, isnr) = n_biterror/(N_packet*L_frame*Nmod);
end
%% 画图
BERs(2,:) = [1.40608250000000	1.27856000000000	1.12480750000000	0.941037500000000	0.738272500000000	0.582960000000000	0.417507500000000	0.287120000000000	0.197902500000000	0.115562500000000	0.0775400000000000];
BERs(3,:) = [0.892967500000000	0.736142500000000	0.576275000000000	0.428812500000000	0.306602500000000	0.212430000000000	0.143212500000000	0.0964100000000000	0.0644250000000000	0.0476550000000000	0.0363025000000000];
BERs(4,:) = [1.16728000000000	1.02494250000000	0.830452500000000	0.641750000000000	0.481362500000000	0.320395000000000	0.211695000000000	0.144177500000000	0.0884050000000000	0.0634100000000000	0.0345100000000000];
BERs(5,:) = [0.679290000000000	0.521507500000000	0.376435000000000	0.253155000000000	0.160840000000000	0.102062500000000	0.0656150000000000	0.0394150000000000	0.0272925000000000	0.0190700000000000	0.0148800000000000];
for i = 1:5
    semilogy(SNRs_dB,BERs(i,:), gss(i));
    hold on;
end