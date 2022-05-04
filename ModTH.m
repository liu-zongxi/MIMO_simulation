%------------------------TH编码中的Mod功能------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年5月4日20点27分-----------------%
function [y]=ModTH(x,A)
% 输入
% x：要被mod的数
% A: mod的基数
% 输出
% y：mod后的值
int_real=floor((real(x)+A)/(2*A));
int_imag=floor((imag(x)+A)/(2*A));
y=x-int_real*(2*A)-1j*int_imag*(2*A);