
function [W, gamma] = init_gama(base_weight,S1,neighbor,lambda,distF)

[~,nBase] = size(S1); %
[nSmp,~]=size(S1{1,1});

W = zeros(nSmp);
gamma_temp = zeros(nSmp,1);
D = zeros(1,nSmp);
k=neighbor;
if k==nSmp
    k=k-2;
end
for i = 1:nSmp
    %按行计算
    row_temp = zeros(1,nSmp);
    for p = 1:nBase
        temp_s=S1{1,p};           
        row_temp = row_temp + base_weight(p) * temp_s(i,:);
    end

    D = (lambda*distF(i,:) - row_temp); % 行向量
%     D =  - row_temp;
    [~, idx] = sort(D, 2); 

    id = idx(1,2:k+2); % 第一个是对角线不取
    di = D(1, id); 
    W(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps); 
    gamma_temp(i) = k/2*di(k+1) - 1/2*sum(di(1:k));
%     gamma_temp(i) = k/2*di(k+1) - 1/2*sum(di(1:k));
end
gamma=gamma_temp;






