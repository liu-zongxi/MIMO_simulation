%-----------------------OSTBC的天线选择---------------------%
%-----------------------author:lzx--------------------------%
%-----------------------date:2022年4月30日09点58分-----------------%
%% 参数设置
clear;
clf;
L_frame = 1000;
Niter = 400;
Nmod = 2;
M = 2^Nmod;

N_T_Alamouti = 4;
NT = 2;
NR = 1;

SNRs_dB = 0:2:20;
SNRs = 10.^(SNRs_dB/10);
N_SNR = length(SNRs);

BERs = zeros(2, N_SNR);
%% 主函数
for isnr = 1:N_SNR
    SNR = SNRs(isnr);
    sigma = sqrt(NT/(2*SNR));
    n_biterror = 0;
    for iiter = 1:Niter
        % 生成数据
        frame_origin = randi([0,1],L_frame,NT*Nmod);
        % QPSK调制
        frame_mod=QPSKMod(frame_origin,L_frame, NT);
        frame_transmit = reshape(frame_mod, NT, 1, L_frame);
        frame_alamouti = [frame_transmit(1,1, :) -conj(frame_transmit(2,1, :));
                          frame_transmit(2,1, :) conj(frame_transmit(1,1, :))];
        % 生成信道，SIMO有NR个信道
        Hiid = (randn(NR,N_T_Alamouti)+1j*randn(NR,N_T_Alamouti))/sqrt(2);
        noise = sigma*(randn(NR,NT)+1j*randn(NR,NT));
        for ih = 1:N_T_Alamouti
            fro_h(ih) = norm(Hiid(:,ih),'fro');
        end
        [val,index] = sort(fro_h,'descend');
        H_sel = Hiid(:,index([1 2]));
        for iframe = 1:L_frame
            y(:,:,iframe) = H_sel*frame_alamouti(:,:,iframe)+noise;
        end
        % 接收
        H_sel_re = [conj(H_sel(1)) H_sel(2); conj(H_sel(2)) -H_sel(1)];
        y_re = [y(1,1,:); conj(y(1,2,:))];
        for iframe = 1:L_frame
            frame_re(:,:,iframe) = (H_sel_re* y_re(:,:,iframe))./(norm(H_sel_re,"fro")^2);
        end
        frame_pre_demod = reshape(frame_re, L_frame, NT);
        frame_demod = QPSKDemod(frame_pre_demod,L_frame,NT);
        % 计算误码率
        n_biterror_tmp = sum(sum(abs(frame_demod - frame_origin)))
        n_biterror = n_biterror + n_biterror_tmp;
    end
    BERs(isnr) = n_biterror/(Niter*L_frame*Nmod*NT);
end
%% 画图
BERs(2, :) = [0.177905000000000	0.131740000000000	0.0911950000000000	0.0580825000000000	0.0335800000000000	0.0178250000000000	0.00892500000000000	0.00408000000000000	0.00184000000000000	0.000750000000000000	0.000295000000000000];
semilogy(SNRs_dB,BERs(1, :),"-^", SNRs_dB,BERs(2, :),"-x");
axis([SNRs_dB([1 end]) 1e-6 1e0]);
xlabel('SNR[dB]'), ylabel('BER'); 
legend('Precoded Alamouti','No Precoded Alamouti');
