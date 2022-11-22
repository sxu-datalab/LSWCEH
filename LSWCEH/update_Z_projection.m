function Z = update_Z_projection(A)
% min_{Z >=0, Z = Z', Z*1=1, Tr(Z) = k}  ||S-V||_F^2
nSmp = size(A, 1);
Z = A;
maxIter = 10;
for ite = 1:maxIter
  
%     %******************************************
    % Lemma 2, KDD, 2016
    %*******************************************
%     Z = Z - diag(diag(Z));
    Z = (Z + Z') / 2;
    Z = Z + (nSmp + sum(sum(Z)))/(nSmp^2) - sum(Z, 2)/nSmp - sum(Z, 1)/nSmp;
% Z = Z + (nSmp+sum(sum(Z)))/(nSmp^2)*ones(nSmp) - 1/nSmp*repmat(sum(Z),nSmp,1) - 1/nSmp*repmat(sum(Z,2),1,nSmp);
    %*******************************************
    % Lemma 3, KDD, 2016
    %*******************************************
    Z = max(Z, 0);
    Z=min(Z,1);
    %Z = Z - diag(diag(Z));
%     dz = project_simplex(diag(Z), c);
%     dz = EProjSimplex_new(diag(Z), c);
%     Z = Z - diag(diag(Z)) + diag(dz);    
end