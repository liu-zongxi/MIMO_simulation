%-----------------------码本生成-------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年4月30日09点58分-----------------%
function [code_book]=CodebookGenerator(NT, M, L)
%输入
%NT: 天线数
%M: 列向量长度，看书码字长度
%L: 码本长度
%输出
% code_book：码本
if NT == 2 && M == 1 && L ==8
    cloumn_index = [1];
    rotation_vector = [1,0];
elseif NT == 3 && M == 1 && L ==32
    cloumn_index=[1]; 
    rotation_vector=[1 26 28];
elseif NT == 4 && M == 2 && L ==32
    cloumn_index=[1 2]; 
    rotation_vector=[1 26 28];
elseif NT == 4 && M ==1 && L == 64
    cloumn_index=[1]; 
    rotation_vector=[1 8 61 45];
elseif NT == 4 && M ==2 && L == 64
    cloumn_index=[0 1]; 
    rotation_vector=[1 7 52 56];
elseif NT == 4 && M ==3 && L == 64
    cloumn_index=[0 2 3]; 
    rotation_vector=[1 8 61 45];
else
    error("没有对应的码本");
end
index_w = 0:NT-1;
w = exp(1j*2*pi/NT*(index_w).'*index_w)/sqrt(NT);
w_1 = w(:, cloumn_index+1);
theta = diag(exp(1j*2*pi/L*rotation_vector));
code_book(:,:,1) = w_1;
for i = 1:L-1
    code_book(:,:,i+1) = theta*code_book(:,:,i);
end