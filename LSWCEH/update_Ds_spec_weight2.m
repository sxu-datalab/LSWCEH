function [S_new,F] = update_Ds_spec_weight2(S,lambda,nCluster)
%UPDATE_DS 此处显示有关此函数的摘要
%   此处显示详细说明
% **************************************************************************
% (1/wi)||S^-Si||^2 + lambda*2*trace(F'*L*F)
% **************************************************************************
% init F
% lambda=1;
converges = 0;
iter=1;
myeps=1e-3;
maxIter=20;

%init S0
nweight=length(S);
beta  = ones(nweight,1)/nweight;
[nSmp,~]=size(S{1,1});
S0=zeros(nSmp,nSmp);
C=sum(1./beta);
for i=1:nweight
    S0=S0+1/beta(i)*S{1,i};
end
% S0=S0./C;
% S0=mapminmax(S0,0,1);
%init F
Das=diag(sum(S0));
L = Das - (S0 + S0')/2;
L = (L + L')/2;
[eigvec, eigval] = eig(L);
[~, eigIdx] = sort(diag(eigval), 'ascend');
F = eigvec(:, eigIdx(1:nCluster));
S_new=S0;

% lambda=nSmp;
objhistory=[];
while ~converges
    S0=zeros(nSmp,nSmp);
    for i=1:nweight
        S0=S0+1/beta(i)*S{1,i};
    end
%     S0=mapminmax(S0,0,1);
    %update S
    distF = EuDist2(F, F, 0);
    lambda=(norm(S0./C,'fro'))/(norm(distF,'fro'));
    B=(S0-lambda*distF)./C;
    B=max(0,B);
%     S_old=S_new;
    S_new=update_Z_projection_32(B, nCluster);

    %updata beta
    e=zeros(nweight,1);
    
    for i1=1:nweight
        d=S_new-S{i1};
        e(i1,1)=norm(d,'fro');
    end   
    beta=  e./sum(e);  
    C=sum(1./beta);
    %update F
    F_old = F;
    Das=diag(sum(S_new));
    L = Das - (S_new + S_new')/2;
    L = (L + L')/2;
    [eigvec, eigval] = eig(L);
    [eigval2, eigIdx] = sort(diag(eigval), 'ascend');
    F = eigvec(:, eigIdx(1:nCluster));
   
%     fn1 = sum(eigval2(1:nCluster));
%     fn2 = sum(eigval2(1:nCluster + 1));
%     if fn1 > 0.00000000001
%         lambda = 2*lambda;
%     elseif fn2 < 0.00000000001
%         lambda = lambda/2;
%         F = F_old;
%     end
    obj=culobj(S_new,S0,lambda,F,L);
    objhistory(iter)=obj; %#ok
    
    if (iter > 2 && abs((objhistory(iter))-(objhistory(iter-1))) /(objhistory(iter-1)) < myeps) || iter > maxIter
        converges = 1;
    end
    iter=iter+1;
%     F = bsxfun(@rdivide, F, sqrt(sum(F.^2, 2)) + eps);
end
end

function obj=culobj(S_new,S_0,lambda,F,L)
o1=lambda*trace(F'*L*F);
o2=norm(S_new-S_0);
obj=o1+2*o2;
end