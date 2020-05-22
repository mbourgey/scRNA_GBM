% roadmap.m
% Creator: Charles Couturier
% Created: January 2018
% Last update: April 2020
% This script loads a dataset containing both fetal and cancer cells, builds a
%  developmental roadmap from selected fetal cell types and overlays the cancer cells 
%  onto it. Based on their position on the roadmap, cell types are then attributed to the 
%  cancer cells. The dpt toolbox needs to be in the path.

rng(1);
post_dim = 10; % number of PCA dimensions to include in analysis

%% load data, calculate cc score, define subgroups
% gbm_ds.mat contains both cancer and fetal cells
load('gbm_ds.mat') % downsampled data with equal cells per sample to avoid bias
% load fetal clustering solution
load('community_final.mat')  
X = matrix;

[G1S,G2M] = get_ccscore(X,genes(:,2),1); % based on Tirosh et al.

cc_thresh = 0;
not_cycling = (G2M<cc_thresh & G1S<cc_thresh);

fetal = 1; % fetal id within 'sample'
% clusters to include in 'clusters': tRG, neuron, astro, GPC, OLC
include_fetal_cluster = [1,6,8,9,10];
include_fetal = logical(zeros(size(X,1),1));
which_fetal = logical(zeros(sum(sample==fetal),1));

% select k cells from each fetal cluster to avoid bias in the roadmap
for i=1:length(include_fetal_cluster)
    Ncells(i) = sum(cic==include_fetal_cluster(i));
end
k = min(Ncells);

for i=1:length(include_fetal_cluster)
    fetali = ( cic==include_fetal_cluster(i) );
    Ni = sum(fetali);
    subsample = zeros(Ni,1);
    subsample(randsample(Ni,k)) = 1;
    which_fetal(fetali) = subsample;
end
include_fetal(sample==fetal) = which_fetal;


% redefine some group names
cancer_id = sample_id(1:numel(sample_id) ~= fetal);
cancer_sample = sample(sample~=fetal);

%fetal_sample = cic(ismember(cic,include_fetal_cluster));
fetal_sample = cic(which_fetal);

% decide which cells to include in PCA
in_cells = (include_fetal & not_cycling');
Xin = X(in_cells,:);

X_fetal = X(include_fetal,:); %sample==fetal
X_cancer = X(sample~=fetal,:);

%% create roadmap. PCA and diffusion are different representations of the same entity.
% PCA
[coeff,score,latent] = pca(Xin);
Xpost = X * coeff(:,1:post_dim);

Xpost = X * coeff(:,1:post_dim);
Xpost_fetal = X_fetal * coeff(:,1:post_dim);
Xpost_cancer = X_cancer * coeff(:,1:post_dim);

% diffusion
k = 50;
nsig = 10;
l = 10;
data_fetal = Xpost_fetal(:,1:l);
data_cancer = Xpost_cancer(:,1:l);
[Tf,Tc,phi0] = diffusionmap.T_nn(data_fetal,data_cancer,k,nsig);
[ phi, lambda ] = diffusionmap.eig_decompose_normalized(Tf,l);
phi_cancer = Tc*phi;

% control
Xpost_control = X_cancer(:,randperm(size(X_cancer,2))) * coeff(:,1:post_dim);
data_control = Xpost_control(:,1:l);
[Tf2,Tcontrol,phi0] = diffusionmap.T_nn(data_fetal,data_control,k,nsig);
[ phi2, lambda2 ] = diffusionmap.eig_decompose_normalized(Tf2,l);
phi_control = Tcontrol*phi2;

%% classify cancer cells into cell types based on roadmap position
id_fetal = knnsearch(Xpost_fetal, Xpost_cancer, 'K', 5);
knn_class = fetal_sample(id_fetal);
cancer_class = mode(knn_class,2);

MdlLinear = fitcdiscr(phi, fetal_sample);
[opti_class,score,~] = predict(MdlLinear,phi_cancer);

cancer_class = zeros(size(opti_class));
%
% top percent
[row, col]=find(score>0.9999); % 0.999
cancer_class(row)=col;


%% parameters for general figures and hierarchy

dim1 = 2;
dim2 = 3;
dim3 = 5;


% KS for control vs cancer distribution in diff
for i = 1:3
    d = i+1;
    [~,KSp(i)] = kstest2(phi_control(:,d),phi_cancer(:,d));
end
fprintf('KS pvalue: %.2e\n',KSp);


%% cell cycle
cycling = (G1S(sample~=fetal)>1.5 | G2M(sample~=fetal)>1.5);
%cycling = (G1S(include_fetal)>1.5 | G2M(include_fetal)>1.5);
data = phi_cancer;
%data = phi;

% figure specs
xpos = 100;
ypos = 50;
width = 800;
height = 1600;
az = 60;
el = 45;

f1=figure;
set(f1, 'Position', [xpos ypos width height])
cmap = lines(4);
types = unique(fetal_sample);
for i = 0:1
    c = i+1;
    scatter3(phi_cancer(cycling==i,dim1),phi_cancer(cycling==i,dim2),-phi_cancer(cycling==i,dim3),15,cmap(c,:),'filled')
    hold on
end 
view(az,el);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'ZTickLabel',[]);

figure;
gscatter(data(:,3),data(:,4),cycling,lines)
xlabel('DC1');ylabel('DC2');
legend({'non-cycling','cycling'})

figure;
S = -data(:,dim3);
nbins = 15;
[N,edges] = histcounts(S(cycling),nbins);
[N0,~] = histcounts(S(~cycling),edges);
centers = (edges(1:end-1) + edges(2:end))/2;
bar(centers, N./(N+N0));
ylim([0 1])
xlabel('GPC score')
ylabel('Proportion of cycling cells')

for i = 1:length(cancer_id)
    x = i+1;
    ccP(i) = sum(cycling(cancer_sample==x))/size(cycling(cancer_sample==x),2);
end

%% heatmaps + hierarchy figure

% figure specs
xpos = 100;
ypos = 50;
width = 800;
height = 1600;
az = 60;
el = 45;

cmap_dum = lines(5);

% fetal
f1=figure;
set(f1, 'Position', [xpos ypos width height])
cmap = [cmap_dum(5,:); cmap_dum(1:4,:)];
types = unique(fetal_sample);
for i = 1:length(types)
    c = i;
    cl = types(i);
    scatter3(phi(fetal_sample==cl,dim1),phi(fetal_sample==cl,dim2),-phi(fetal_sample==cl,dim3),15,cmap(c,:),'filled')
    hold on
end  
view(az,el);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'ZTickLabel',[]);

% cancer
f2=figure;
set(f2, 'Position', [xpos ypos width height])
g = 0.5;
cmap2 = [g g g; cmap];
for i = 0:max(cancer_class)
    c = i+1;
    scatter3(phi_cancer(cancer_class==i,dim1),phi_cancer(cancer_class==i,dim2),-phi_cancer(cancer_class==i,dim3),10,cmap2(c,:),'filled')
    hold on
end    
view(az,el);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'ZTickLabel',[]);

% cancer by patient
figure;
cancer_prog = cancer_class==4;
cmap_prog = cmap2([1 5],:);
for i=2:13
    subplot(4,3,i-1)
    for j = [0 1]
        c = j+1;
        plot_cell = cancer_sample==i & cancer_prog==j;
        scatter3(phi_cancer(plot_cell,dim1),phi_cancer(plot_cell,dim2),-phi_cancer(plot_cell,dim3),10,cmap_prog(c,:),'filled')
        hold on
    end 
    view(az,el);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'ZTickLabel',[]);
end

% cancer heatmap with label
%top = 1000; %500;
ngenes = 200; %500
[rho,pval] = corr(phi_cancer(cancer_class~=0,:),clogm(cancer_class~=0,:)); 

Tgenes = table();
for i = [dim1 dim2 dim3] %size(phi_cancer,2)
    [~,gI] = sort(rho(i,:),'descend');
    [~,dI] = sort(phi_cancer(:,i),'ascend');
    geneI = [gI(1:ngenes) gI(end-ngenes+1:end)]; 
    data2 = cmatrix(dI,geneI);
    
    HeatMap(data2','Colormap',redbluecmap,'DisplayRange',3);
    %dumT = table(cgenes(gI,2),pval(i,gI)');
    %dumT.Properties.VariableNames = {sprintf('DC%i',i), sprintf('p%i',i)};
    dumT = table(cgenes(gI,2));
    dumT.Properties.VariableNames = {sprintf('DC%i',i)};
    Tgenes = [Tgenes dumT];
    
    figure;
    class_plot = cancer_class(dI)';
    
    imagesc(class_plot)
    colormap(cmap2)
    axis off
end
writetable(Tgenes,'gene_order_by_DC.csv')

    
%% assess TCGA
diff2D = [1,2];
point_size = 5;
include_cell_type = [2,3,4,5,6];

TCGA = categorical(TCGA_class);
TCGA = categorical(TCGA, categories(TCGA), TCGA_subtype);

celltype = categorical(cancer_class);
nsubtype = [{'unspecified'}, clusters(include_fetal_cluster)];
celltype = categorical(celltype, categories(celltype), nsubtype);

for j=1:length(cancer_id)
    s = j+1; % add one for correct sample
    tcga = TCGA(cancer_sample==s);
    ct = celltype(cancer_sample==s);
    for i=1:length(include_cell_type)
        cTCGA(:,i,j) = histcounts( tcga(ct==nsubtype(include_cell_type(i))) );
    end

    for i=1:length(TCGA_subtype)
        cCT(:,i,j) = histcounts( ct(tcga'==TCGA_subtype(i)) );
    end
end
cCT = cCT(include_cell_type,:,:);

figure;
x0=500;
y0=500;
width=600;
height=600;
set(gcf,'units','points','position',[x0,y0,width,height])
plot_rows = ceil(length(include_cell_type)/2);
for i=1:length(include_cell_type)
    subplot(plot_rows,2,i)
    pie(sum(cTCGA(:,i,:),3))%,TCGA_subtype)
    title(sprintf('%s',nsubtype{include_cell_type(i)}))
end


%% combined progenitor score

% fetal
Sf = Xpost_fetal(:,2)-Xpost_fetal(:,4);
Lf = Xpost_fetal(:,1);
figure;
gscatter(Lf,Sf,fetal_sample,cmap)
xlabel('Lineage');ylabel('Progenitor');
legend(clusters(include_fetal_cluster))

% cancer
S = Xpost_cancer(:,2)-Xpost_cancer(:,4);
L = Xpost_cancer(:,1);
figure;
scatter(L,S,5,'filled')
xlabel('Lineage');ylabel('Progenitor');

