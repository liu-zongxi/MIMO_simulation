%-----------------------MIMO信道的CDF----------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:09点34分-----------------%
% Ergodic_Capacity_CDF.m
%% 设置参数
clear;
figure
SNR_dB = 10; % 设定信噪比
SNR = 10.^(SNR_dB / 10.); % 转化信噪比
N_iter = 50000; % 迭代次数
grps = ['b:'; 'b-']; % 画图
C = zeros(1, N_iter); % 信道容量初始化
N_hist = 50; % 直方图分成多少份
CDFs = zeros(2, N_hist); % CDF 初始化
Rates = zeros(1, N_hist); %传输速率,这实际上就是容量

%% 主函数
for Icase = 1:2
    % 测试2*2和4*4两种情况
    if Icase == 1
        nT = 2;
        nR = 2; % 2x2
    else
        nT = 4;
        nR = 4; % 4x4
    end

    rank = min(nT, nR); % 秩
    I = eye(rank);

    for iter = 1:N_iter
        H = sqrt(1/2) * (randn(nR, nT) + 1j * randn(nR, nT)); % 先假设信道是完全独立的，信道就可以建模为瑞利信道
        C(iter) = log2(real(det(I + SNR / nT * (H' * H))));% 信道容量计算，H'*H本身计算结果就是实数，这里只是做一个类型转换
        % C(iter) = log2(det(I + SNR / nT * (H' * H)));% 信道容量计算
    end
    figure(1);
    hist = histogram(C, N_hist);
    PDF = hist.Values / N_iter;

    for i = 1:N_hist
        Rates(i) = (hist.BinEdges(i) + hist.BinEdges(i + 1)) / 2;
    end

    for i = 1:N_hist
        CDFs(Icase, i) = sum(PDF([1:i]));
    end

    figure(2);
    plot(Rates, CDFs(Icase, :), grps(Icase, :));
    hold on
end

%% 画图
xlabel('Rate[bps/Hz]');
ylabel('CDF');
axis([1 18 0 1]);
grid on;
set(gca, 'fontsize', 10);
legend('{\it N_T}={\it N_R}=2', '{\it N_T}={\it N_R}=4');
