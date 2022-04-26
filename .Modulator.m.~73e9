%------------------瑞利信道下MRC分集的性能----------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年4月22日23点46分-----------------%
% 这个调制函数非常精彩，决定在抄一下
function [symbol_mod,sym_table,M] = Modulator(input_data, nmod)
N_data = length(input_data);
if nmod==1      % BPSK modulation
    sym_table=exp(1j*[-pi 0]);    %BPSK 0和1对应-1和1
    input_symbol=input_data;     % 只有一位，直接输入
    symbol_mod=sym_table(input_symbol+1);  
    M=2;
elseif nmod==2
    sym_table = exp(1j*pi/4*[-3 3 -1 1]);
    input_symbol = reshape(input_data,nmod,N_data/nmod);
    symbol_mod=sym_table([2 1]*input_symbol+1);    %括号内是二进制运算，加一对应index
    M=4;
elseif nmod==3
    sym_table = exp(1j*pi/4*[0 1 3 2 6 7 5 4]);
    input_symbol = reshape(input_data,nmod,N_data/nmod);
    symbol_mod=sym_table([4 2 1]*input_symbol+1);    %括号内是二进制运算，加一对应index
    M=8;
elseif nmod==4
    m=0;
    sym_table = zeros(1,16);
    for k=-3:2:3
      for l=-3:2:3
         m=m+1;
         sym_table(m) = (k+1j*l)/sqrt(10); % power normalization
      end
    end
    sym_table = exp(1j*pi/4*[0 1 3 2 6 7 5 4]);
    input_symbol = reshape(input_data,nmod,N_data/nmod);
    symbol_mod=sym_table([8 4 2 1]*input_symbol+1);    %括号内是二进制运算，加一对应index
    M=16;
else
    error('Unimplemented modulation');
end
