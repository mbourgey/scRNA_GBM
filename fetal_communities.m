% fetal_communities.m
% Creator: Charles Couturier
% Date created: January 2018
% Last update: April 2020
% This script clusters scRNA cells using community detection.
% BCT and dpt toolboxes need to be in path.

import diffusionmap.*;
import dpt.*;
import third_party.*;
% run community detection on bsr matrix
load('gbm_visualize.mat') % a file containing the data as well as its PCA and tSNE projections

%% edge matrix

data = Xpost; % Xpost contains the data in the form of its top 10 PCA components
k = 50;

% knn euclidean distance
ch_s=min((size(data,1)-1),1000); % chunk size is 1000
lnn = parfor_spdists_knngraph(data, k, 'distance', 'Euclidean', 'chunk_size', ch_s, 'verbose', true);
lnn = spdists_undirected(lnn); %make lnn undirected
cij = 1./(1+lnn);

% dpt, for a less sparse visualization of cell proximity later

k2 = k;
nsig = 10;
[T, ~,phi0] = T_nn(data,data,k2,nsig); %[T, ~,phi0] = T_loc(data,data,nsig);
n = size(T, 1);
M = (eye(n) - T + phi0 * phi0' )^-1 - eye(n); % accumulated transition matrix
D = squareform(pdist(M));
cij2 = 1./(1+D);


% keep knn only 
%cij( (lnn+eye(size(lnn)))==0 ) = 0;  %exclude 0 values for lnn except self
cij( (lnn)==0 ) = 0;

clearvars T D M

%% general parameters
n = size(cij,1);

gamvals = 0:0.01:2;%0.8:0.02:1.3; % gamma values over which to iterate
ngam = length(gamvals);

nreps = 100;
ci = zeros(n,ngam,nreps);   % community assignments
q = zeros(ngam,nreps);      % modularity values
ciu = zeros(n,ngam);        % consensus communities

parpool; % begin parallelization

%% community detection
for igam = 1:ngam
    gam = gamvals(igam);
    
    fprintf('gam %1.2f ',gam);
    tic;
    parfor irep = 1:nreps
        [ci(:,igam,irep),q(igam,irep)] = community_louvain(cij,gam,[],[]); % 'negative_sym' or 'negative_asym' or []
        %[ci(:,igam,irep),q(igam,irep)] = modularity_louvain_und_sign(B,'sta');
    end
    
    fprintf(' finished in %.2f s\n',toc);
end

% look at partition similarity
rzvals = zeros(ngam,2);
for igam = 1:ngam
    rz = zeros(nreps);
    mask = triu(ones(nreps),1) > 0;
    for i = 1:(nreps - 1)
        parfor j = (i + 1):nreps
            % z-score of Rand index (similarity between partitions)
            rz(i,j) = fcn_randz(ci(:,igam,i),ci(:,igam,j));
        end
    end
    rzvec = rz(mask);
    rzvals(igam,1) = mean(rzvec);
    rzvals(igam,2) = std(rzvec)^2;
    fprintf('gam %i of %i\n',igam,ngam);
end
snr = rzvals(:,1)./sqrt(rzvals(:,2));
mode_nclusters = mode(max(ci,[],1),3);

% get consensus (see Bassett et al., 2013, Chaos)
bestgam = 11; % input best gamma

citemp = squeeze(ci(:,bestgam,:));
ag = agreement(citemp) / nreps; % agreement (probability)
    
cinull = zeros(size(citemp));
% get null agreement
for irep = 1:nreps
    cinull(:,irep) = citemp(randperm(n),irep);
end
agnull = agreement(cinull) / nreps;
tau = max(agnull(:));
    
if max(citemp(:))>1 % consensus only if >1 partition
    cic = consensus_und(ag,tau,10);
else
    cic = ci(:,bestgam,1);
end

%% save

save 'community.mat' gamvals rzvals cij ci cic bestgam 


%% plots

% plot rzvals
figure;
plotyy(gamvals,rzvals(:,1),gamvals,sqrt(rzvals(:,2)))
legend({'mean randz','sd randz'})

% plot snr
figure;
ax = plotyy(gamvals,snr,gamvals,mode_nclusters);
set(ax(1),'Ylim',[0, max(snr)])
legend({'snr','number of clusters'})

% show matrix - this is where diffusion is helpful
figure;
[x,y,~] = grid_communities(cic); % call function
[~,idx] = sortrows([cic mean(abs(cij2))'],[1 -2]);
%spy(cij(idx,idx))
imagesc(cij2(idx,idx));          % plot ordered adjacency matrix
hold on;                        % hold on to overlay community visualization
plot(x,y,'r','linewidth',2);  
colorbar
axis square

% show clusters on t-SNE plot, assumes tSNE has already been computed and stored
%  in mappedX
figure;
psize = 10;
cmap = cbrewer('qual','Paired',12,'PCHIP');
gscatter(mappedX(:,1),mappedX(:,2),cic,cmap,'.',psize)
ylabel('t-SNE1')
xlabel('t-SNE2')


