function [Z] = update_Zs(base_weight,S1,lambda,distF,gamma)

[~,nBase] = size(S1); %
[nSmp,~]=size(S1{1,1});

for i = 1:nSmp
    %按行计算
    row_temp = zeros(1,nSmp);
    for p = 1:nBase
        temp_s=S1{1,p};           
        row_temp = row_temp + base_weight(p) * temp_s(i,:);
    end

    ft = lambda*distF(i,:) - row_temp;
%     ft =  - row_temp;
    Z_hat = -ft/(2*(gamma(i))+eps);
    
    indx = 1:nSmp;
    indx(i) = [];
    [Z(i,indx), ~] = EProjSimplex_new(Z_hat(:,indx));
end

end

