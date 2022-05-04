%-----------------------Alamouti预编码对比-------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年4月30日09点58分-----------------%
%% 参数设置
clear; 
clf;
NT=4;
NR=4;
N_iter=1000;
I=eye(NR,NR);
gss=['-kx';'-ro';'-k^';'-g+'];
SNRs_dB=0:2:20;
SNRs = 10.^(SNRs_dB/10);
N_SNR = length(SNRs);
Capacitys = zeros(1,N_SNR);
%% 主函数
for n_ant=1:NT
    if n_ant>NT || n_ant<1
       error('n_ant must be between 1 and NT!');
     else
        Types = nchoosek([1:NT],n_ant);
        N_type = size(Types,1);
     end
   for isnr=1:N_SNR
      SNR = SNRs(isnr)/n_ant;  
%       rand('seed',1); randn('seed',1);  
      sum_capacity = 0;
      for iiter=1:N_iter
         H = (randn(NR,NT)+1j*randn(NR,NT))/sqrt(2);
         log_SH = zeros(1,N_type);
         for itype=1:N_type
            H_sel = H(:,Types(itype,:)); 
            log_SH(itype)=log2(real(det(I+SNR*H_sel*H_sel'))); % Eq.(12.22)
         end
         sum_capacity = sum_capacity + max(log_SH);
      end
      Capacitys(isnr) = sum_capacity/N_iter;
   end
   plot(SNRs_dB,Capacitys,gss(n_ant,:), 'LineWidth',2); hold on;
end
xlabel('SNR[dB]'), ylabel('bps/Hz'), grid on;
legend('sel-ant=1','sel-ant=2','sel-ant=3','sel-ant=4')
