function [label1,label2,label3,label4] = DAMCEL(baseCls, nCluster,lambda)

[nSmp, nBase] = size(baseCls);
Hs = baseCls2Hs(baseCls);

if ~exist('lambda', 'var')
    lambda = 1;
end

iter=1;
myeps=1e-3;
% **************************************************************************
% S1:A
% **************************************************************************
% S1=cell(1, nBase);
base_weight=ones(nBase,1)/nBase;
G=zeros(nSmp,nCluster);
Q=cell(1,nBase);
for i1 = 1:nBase
    G=G+base_weight(i1,1)*Hs{i1};
    Q{i1}=eye(nCluster,nCluster);
end
converge=0;
H=G;

while  ~converge
    %update w   
    z = zeros(nBase, 1);
    for i1 = 1:nBase  
        z(i1) = trace(H'*Hs{i1}*Q{i1});
    end
    base_weight = z / sqrt(sum(z.^2) + eps);
%     base_weight=e./(norm(e,2));
    G=zeros(nSmp,nCluster);
    for i1 = 1:nBase
    G=G+base_weight(i1,1)*Hs{i1};
    end
    % updata H
    H_old=H;
    distF=EuDist2(H_old, H_old, 0);
    H=update_H(H_old,distF,lambda,nCluster,G);
    %update Q
    for i1 = 1:nBase
        P_i=H'*Hs{i1};
        [U, ~, V] = svd(P_i);
        Q{i1} = U * V';
    end
       
    %update F
    S=H*H';
    S=(S + S')/2;
    Das=diag(sum(S));
    L = Das - (S + S')/2;
    [eigvec, eigval] = eig(L);
    [eigval2, eigIdx] = sort(diag(eigval), 'ascend');
    F = eigvec(:, eigIdx(1:nCluster));
    fn1 = sum(eigval2(1:nCluster));
    fn2 = sum(eigval2(1:nCluster + 1));
    F_old = F;
    if fn1 > 0.00000000001
        lambda = 10*lambda;
    elseif fn2 < 0.00000000001
        lambda = lambda/10;
        F = F_old;
    elseif rank(L)==nCluster
        break
    end
    iter=iter+1;
    if iter>2 && norm((H-H_old),'fro')<myeps || iter>10
        converge=1;
    end
end

S=H*H';
S=(S + S')/2;
label1=litekmeans(H, nCluster, 'MaxIter', 100, 'Replicates', 20);
% label1 = litekmeans(H, nCluster, 'emptyaction', 'singleton', 'replicates', 50, 'display', 'off');
CKSym = BuildAdjacency(S,nCluster);
label2 = SpectralClustering(CKSym,nCluster,3);
Das=diag(sum(S));
L = Das - S;
L = (L + L')/2;
[P, ~] = eigs(L, nCluster, 'smallestreal');
P = bsxfun(@rdivide, P, sqrt(sum(P.^2, 2)) + eps);
% label3 = litekmeans(P, nCluster, 'emptyaction', 'singleton', 'replicates', 50, 'display', 'off');
label3=litekmeans(P, nCluster, 'MaxIter', 100, 'Replicates', 20);
[clusternum, label4] = graphconncomp(sparse(S));
label4 = label4';
if clusternum ~= nCluster
    sprintf('Can not find the correct cluster number: %d', nCluster);
end
end

function Hs = baseCls2Hs(baseCls)
[nSmp, nBase] = size(baseCls);
Hs = cell(1, nBase);
for iBase = 1:nBase
    [~, ~, label] = unique(baseCls(:, iBase));
    nc = max(label);
    Xi = zeros(nSmp, nc);
    lidx = sub2ind([nSmp, nc], (1:nSmp)', label);
    Xi(lidx) = 1;
    Hs{iBase} = Xi;
end
end
