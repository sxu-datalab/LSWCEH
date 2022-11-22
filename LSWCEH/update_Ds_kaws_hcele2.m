function [S_new,F] = update_Ds_kaws_hcele2(S,lambda,nCluster)
%UPDATE_DS 此处显示有关此函数的摘要
%   此处显示详细说明
% **************************************************************************
% wi||S^-Si||^2 + lambda*2*trace(F'*L*F) 
% **************************************************************************
% init F
% lambda=1;
converges = 0;
iter=1;
myeps=1e-3;
maxIter=10;
%init S0
nweight=length(S);
beta  = ones(nweight,1)/nweight;
[nSmp,~]=size(S{1,1});
S0=zeros(nSmp,nSmp);
for i=1:nweight
    S0=S0+beta(i)*S{1,i};
end
%init F
Das=diag(sum(S0));
L = Das - (S0 + S0')/2;
L = (L + L')/2;
[eigvec, eigval] = eig(L);
[~, eigIdx] = sort(diag(eigval), 'ascend');
F = eigvec(:, eigIdx(1:nCluster));
S_new=S0;
delta=10;
while ~converges
    S0=zeros(nSmp,nSmp);
    for i=1:nweight
        S0=S0+beta(i)*S{1,i};
    end
    %update S
    distF = EuDist2(F, F, 0);
    B=S0-lambda*distF;
    S_old=S_new;
    S_new=update_Z_projection(B);
    %updata beta
    e=zeros(nweight,1);
    for i1=1:nweight
        e(i1,1)=norm(S{i1}- S_new,'fro')^2;
    end
    e_ba=mean(e);
    beta=exp((-delta*e./e_ba))./sum(exp((-delta*e./e_ba)));

    %update F
    F_old = F;
    Das=diag(sum(S_new));
    L = Das - (S_new + S_new')/2;
    L = (L + L')/2;
    [eigvec, eigval] = eig(L);
    [eigval2, eigIdx] = sort(diag(eigval), 'ascend');
    F = eigvec(:, eigIdx(1:nCluster));
    fn1 = sum(eigval2(1:nCluster));
    fn2 = sum(eigval2(1:nCluster + 1));
    if fn1 > 0.00000000001
        lambda = 2*lambda;
    elseif fn2 < 0.00000000001
        lambda = lambda/2;
        F = F_old;
    end
    
    iter=iter+1;
    if (iter > 10 && ((norm(S_new - S_old) /norm(S_old)) < myeps)) || iter > maxIter
        converges = 1;
    end
    
end
end

