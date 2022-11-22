function [label1,label2] = LSWCEH(baseCls, nCluster,lambda,k,kz,tz)

[nSmp, nBase] = size(baseCls);
Hs = baseCls2Hs(baseCls);


if ~exist('lambda', 'var')
    lambda = 1;
end
if k==1
    neighbor=ceil(nSmp/nCluster);
else
    neighbor=k;
end
% if neighbor>50
%     neighbor=50;
% end
iter=1;
myeps=1e-3;
% **************************************************************************
% S1:A
% **************************************************************************
S1=cell(1, nBase);
for i1 = 1:nBase
    temp=Hs{i1} * Hs{i1}';
%     temp=temp./sqrt(sum(temp,2));
    S1{1,i1} = temp;
    %     S1{1,i1}=S1{1,i1}./sum(S1{1,i1});
end

% **************************************************************************
% init Z gama F
% **************************************************************************

base_weight=sqrt(ones(nBase,1)/nBase);
converge=0;
%init F
init_distF=zeros(nSmp,nSmp);
for i1=1:nBase
    init_distF=init_distF+base_weight(i1)*S1{1,i1};
end

Das=diag(sum(init_distF));
L = Das - (init_distF + init_distF')/2;
L = (L + L')/2;
[F, ~] = eigs(L,nCluster,'smallestreal');
 distF=EuDist2(F, F, 0);


%init gama and Z
[Z,gamma]=init_gama(base_weight,S1,neighbor,lambda,distF);
% gamma=sqrt(ones(nSmp,1).*1e-4);
while  ~converge
    %update w   
    e=zeros(1,nBase);
    for i1=1:nBase
        e(1,i1)=trace(S1{1,i1}*Z);
    end
    base_weight=e./(norm(e,2));
    %update Z
    distF=EuDist2(F, F, 0);
    Z_old=Z;
    Z=update_Zs(base_weight,S1,lambda,distF,gamma);
   
    %update F

    Das=diag(sum(Z));
    L = Das - (Z + Z')/2;
    L = (L + L')/2;
%     [F, ~] = eigs(L,nCluster,'smallestreal');
    [eigvec, eigval] = eig(L);
    [eigval2, eigIdx] = sort(diag(eigval), 'ascend');
    F = eigvec(:, eigIdx(1:nCluster));
    fn1 = sum(eigval2(1:nCluster));
    fn2 = sum(eigval2(1:nCluster + 1));
    F_old = F;
    if fn1 > 0.00000000001
        lambda = 2*lambda;
    elseif fn2 < 0.00000000001
        lambda = lambda/2;
        F = F_old;
    end
    iter=iter+1;
    if iter>2 && norm((Z-Z_old),'fro')<myeps || iter>10
        converge=1;
    end
end
Z=SGD_heat(Z,tz,5);
Z=Z{1,1};
Das=diag(sum(Z));
L = Das - (Z + Z')/2;
L = (L + L')/2;
[P, ~] = eigs(L, nCluster, 'smallestreal');
P = bsxfun(@rdivide, P, sqrt(sum(P.^2, 2)) + eps);
label1 = kmeans(P, nCluster, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off');
CKSym = BuildAdjacency(Z,nCluster);
label2 = SpectralClustering(CKSym,nCluster,3);
% label2=litekmeans(P,nCluster, 'MaxIter',100, 'Replicates',100);

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
