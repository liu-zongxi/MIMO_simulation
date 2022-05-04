%-----------------------降低复杂度的天线选择-------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年4月30日09点58分-----------------%
%% 参数设置
clear; clf;
method_sel = 1; % 0/1 for increasingly/decreasingly ordered selection

NT = 5;
NR = 4;
I = eye(NR,NR);
SNRs_dB = 0:2:20;
SNRs = 10.^(SNRs_dB/10);
N_SNR = length(SNRs);
Niter = 1000;
ns_antenna_select= [1 2 3 4];
N_ant = length(ns_antenna_select);
gss=['-ko';'-k^';'-kd';'-ks'];

%% 主函数
for iant = 1:N_ant
    for isnr = 1:N_SNR
        SNR = SNRs(isnr)/ns_antenna_select(iant);
        sum_capacity = 0;
        for i = 1:Niter
            if method_sel==0
                sel_ant_indices=[];
                rem_ant_indices=1:NT;
            else
                sel_ant_indices=1:NT;
            end
            H = (randn(NR,NT)+1j*randn(NR,NT))/sqrt(2);
            if method_sel == 0
                for n_antenna_select = 1:ns_antenna_select(iant)
                    clear log_SH;
                    for ii = 1:length(rem_ant_indices)
                        H_sel = H(:,[sel_ant_indices rem_ant_indices(ii)]);
                        log_SH(ii) = log2(real(det(I+SNR*H_sel*H_sel')));
                    end
                    [max_cap, max_index] = max(log_SH);
                    sel_ant_index = rem_ant_indices(max_index);
                    rem_ant_indices = [rem_ant_indices(1:max_index-1) rem_ant_indices(max_index+1:end)];
                    sel_ant_indices = [sel_ant_indices sel_ant_index];
                end
            else
                for n_antenna_select = 1:NT-ns_antenna_select(iant)
                    clear log_SH;
                    for ii=1:length(sel_ant_indices)
                        H_sel = H(:,[sel_ant_indices(1:ii-1) sel_ant_indices(ii+1:end)]);
                        log_SH(ii) = log2(real(det(I+SNR*H_sel*H_sel')));
                    end
                    [max_cap, max_index] = max(log_SH);
                    sel_ant_indices = [sel_ant_indices(1:max_index-1) sel_ant_indices(max_index+1:end)];
                end
            end
            sum_capacity = sum_capacity + max_cap;
        end
        sel_capacity(iant,isnr) = sum_capacity/Niter;
    end
    plot(SNRs_dB,sel_capacity(iant,:),gss(ns_antenna_select(iant),:), 'LineWidth',2); 
    hold on;
end