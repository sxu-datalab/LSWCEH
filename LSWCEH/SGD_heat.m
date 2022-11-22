function Ss = SGD_heat(Gss, t, k)
% [1] Diffusion Improves Graph Learning, NIPS, 2019
% [2] Learning with Local and Global Consistency, NIPS, 2003
% [3] The heat kernel as the pagerank of a graph-PNAS-2007
if ~exist('t', 'var')
    t = 1;
end
if ~exist('k', 'var')
    k = 5;
end
Gs=cell(1);
Gs{1}=Gss;
nGraph = length(Gs);
nSmp = size(Gs{1}, 1);
Ss = cell(1, nGraph);
for i1 = 1:nGraph
    %*********************************************
    % Step 1: The symmetric transition matrix
    %*********************************************
    Gsym = (Gs{i1} + Gs{i1}')/2;
    DGsym = 1./sqrt(max(sum(Gsym, 2), eps));
    Gnorm = (DGsym * DGsym') .* Gsym;
    Gnorm = (Gnorm + Gnorm')/2;
    
    %*********************************************
    % Step 2: Heat Kernel diffusion with close-form
    %*********************************************
    S = expm(- t * (eye(nSmp) - Gnorm));
%     S = S - 1e8 * eye(nSmp);
    [Val, Idx] = sort(S, 2, 'descend');

    Idx = Idx(:, 1:k);
    Val = Val(:, 1:k);
    SS = zeros(nSmp);
    
    for iSmp = 1:nSmp
        idxa0 = Idx(iSmp, :);
        ad = Val(iSmp, :);
        SS(iSmp, idxa0) = EProjSimplex_new(ad);
    end
    Ss{i1} = SS;
end
end