% cNMF_seperate.m
% use NMR to find signatures in scRNA data sample by sample. Based on cNMF (Kotliar et al.)
% Author: Charles Couturier
% Date: August 3, 2018

%% prepare data and initialize parameters
% genes: list of gene names and ensembl code; logm: log of raw counts; sample: vector 
%  identifying cells by sample; sample_id: a list of the samples
load('gbm.mat','genes','logm','sample','sample_id')
kvals = 2:15;
L = length(kvals);
nreps = 100;

allHs = cell(L,1);
allHc = cell(L,1);
allSil = cell(L,1);
allE = cell(L,1);
allR = cell(L,1);
allCL = cell(L,1);

for s = 1:length(sample_id)
fprintf('Sample %s in progress, %i of %i\n',sample_id{s},s,length(sample_id))
    
data = logm(sample==s,:);

% initialize
m = size(data,1);
n = size(data,2);

E = zeros(L,1);
R = zeros(L,1);
Sil = zeros(L,1);
consensus_H = cell(L,1);
allH = cell(L,1);
CL = cell(L,1);

%% NMF
opt = statset('MaxIter',1e6);
parfor kid = 1:L
    k = kvals(kid);
    hi = zeros(k*nreps,n);
    consensus_hi = zeros(k,n);
    D = zeros(nreps,1);
    
    for i = 1:nreps
        [~,hi(i*k-k+1:i*k,:),D(i)]=nnmf(data,k,'replicates',20,'options',opt,'algorithm','mult');
    end
    
    % k clusters
    cl = uniform_kmeans(hi,k,'Replicates',10); 
    % find replicates with 1 of each cluster
    goodrep = all(sort(reshape(cl,k,nreps),1) == [1:k]',1); 
    goodrep = repmat(goodrep,k,1);
    goodrep = goodrep(:);
    
    % consensus hi finds median of each component cluster, removing
    % bad replicates (see above)
    for i = 1:k
        consensus_hi(i,:) = median(hi(cl==i & goodrep,:),1); 
    end
    
    allH{kid} = hi;
    consensus_H{kid} = consensus_hi;
    Sil(kid) = mean(silhouette(hi,cl));
    E(kid) = mean(D);
    R(kid) = sum(goodrep)/(k*nreps);
    CL{kid} = cl;
    fprintf('Run for %i components completed\n',k)
end

allHs{s} = allH;
allHc{s} = consensus_H;
allSil{s} = Sil;
allE{s} = E;
allR{s} = R;
allCL{s} = CL;

end

save 'NMF_results_separate.mat' allHs allHc allSil allE allR allCL

%% find best ks

best_k = zeros(length(sample_id),1);
th = 0.9; % threshold in reproducibility to reach to select k

for s = 1:length(sample_id)
    kid = find(kvals==best_k(s));
    E = allE{s};
    R = allR{s};
    Sil = allSil{s};
    CL = allCL{s};
    allH = allHs{s};

    kid = max(find(R>=th));
    best_k(s) = kvals(kid);
    
    figure;
    subplot(1,2,1)
    %subplot(length(sample_id),2,2*s-1)
    yyaxis left
    plot(kvals,E)
    yyaxis right
    plot(kvals,R)
    hold on
    plot(kvals,Sil)
    %legend({'E','Reproducibility of distinct clusters','Silhouette score'})
    
    y = tsne(allH{kid});
    subplot(1,2,2)
    gscatter(y(:,1),y(:,2),CL{kid},lines)
    legend('off')
    title(sprintf('%s',sample_id{s}))
    
    print('-depsc2',sprintf('%s.eps',sample_id{s}))
end


%% get W from H

H = [];
sig_id = [];
for s = 1:length(sample_id)
    
    data = logm(sample==s,:);
    consensus_H = allHc{s};
    kid = find(kvals==best_k(s));
    H0 = consensus_H{kid};
    W0 = max(0,data/H0);

    opt = statset('Maxiter',1e6,'Display','final');
    [~,H1] = nnmf(data,best_k(s),'h0',H0,'w0',W0,...
                     'options',opt,'replicates',100,...
                     'algorithm','als');
    H = [H; H1];
    sig_id = [sig_id; ones(size(H1,1),1)*s];
             
end

save all_signatures.mat H sig_id

