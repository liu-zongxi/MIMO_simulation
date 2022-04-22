%-----------------------频域注水法-------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:09点34分-----------------%
function [gamma] = WaterFilling(H, rank, SNR, nT)
% 输入
% H：MIMO信道
% rank: 信道秩的大小
% SNR：信噪比
% nT:发射天线个数
% 输出
% gamma，注水法生成的系数

sigma = svd(H'*H);      % H'*H的奇异值，因为是艾爾特弥矩阵，奇异值和特征值相等
gamma = zeros(1, rank); % 注水法生成的gamma
index = 1:rank;         % 使用的天线的编号,初始化为全都使用
p=1;
while  p < rank
    index_used = [1:rank-p+1].';    % 被使用的天线的个数
    temp= sum(1./sigma(index(index_used)));
    mu = nT/(rank-p+1.)*(1+1/SNR*temp); % 计算mu
    gamma(index(index_used)) = mu-nT./(SNR*sigma(index(index_used)));   % 计算gamma
    if min(gamma(index))<0     %如果有＜0的结果，这个天线信道应该放弃使用，而把功率重新分配
        i=find(gamma==min(gamma)); % 找到＜0的index
        ii=find(index==i);          % 去除这个天线
        index_new=[index([1:ii-1]) index([ii+1:end])];
        clear index;
        index=index_new;
        p=p+1;
        clear gamma;
        gamma = zeros(1, rank);
    else
        p=rank;                    % 没有的时候就结束
    end
end