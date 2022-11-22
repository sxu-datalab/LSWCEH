function H1=update_H(H_old,D,lambda,nCluster,G)
mu=1;
p=2;
maxiter=10;
converge=0;
H2=H_old';
iter=1;
objhistory=[];
while  ~converge
    B=(H2'-(lambda/mu)*D*H2');
    [U, ~, V] = svd(B);
    H1 = U(:, 1:nCluster) * V';
    H1=max(H1,0);
%     H1=min(H1,1);

    C=(-G+2*lambda*D*H1-2*mu*H1);
    [U, ~, V] = svd(-C);
    H2 = (U(:, 1:nCluster) * V')';

    mu=p*mu;
    iter=iter+1;
    if iter>maxiter
        converge=1;
    end
end
end
