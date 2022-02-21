%The function performs the Differential Expression Analysis of the
%'l2_neurons_olfactory.loom' data published in Zeisel et al (Zeiselet al., 2018)

%(C) Anzhelika Koldaeva, 2021


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%l2_neurons_olfactory.loom   
data = h5read('l2_neurons_olfactory.loom','/matrix');
gene_names = h5read('l2_neurons_olfactory.loom','/row_attrs/Gene');
barcodes = h5read('l2_neurons_olfactory.loom','/col_attrs/CellID');

for i =1:length(barcodes)
    varNames{i} = ['Cell' num2str(i)];
end

data = log2(data+1);
count_matr = array2table(data');
count_matr.(width(count_matr)+1) = gene_names;
count_matr = [count_matr(:,end) count_matr(:,1:end-1)];
count_matr.Properties.VariableNames = ['GeneName' varNames];

%convert gene names
for i = 1:length(gene_names)
    count_matr.GeneName{i} = count_matr.GeneName{i}(find(count_matr.GeneName{i}));
end

genes_to_check = {'Tbx21','Slc32a1','Slc17a7', 'Slc17a6', 'Pkib', 'Cck'}; 
for i = 1:length(genes_to_check)
idx_genes_to_check(:,i) = ismember(count_matr.GeneName, genes_to_check(i));
matr_genes_to_check(i,:) = count_matr{idx_genes_to_check(:,i),2:end};
end

HMobj = clustergram(matr_genes_to_check,... %'ColumnLabels',top500overd_Genes.(1)',...
        'RowLabels',genes_to_check,...
        'Symmetric',false,'Colormap',redbluecmap,'DisplayRange',max(max(matr_genes_to_check)));

%filter genes with respect to the marker genes
indTbxSlc(1) = 1;
for i = 2:size(count_matr,2)
    if count_matr{ismember(count_matr.GeneName, 'Tbx21'),i} > 0 && ...
       count_matr{ismember(count_matr.GeneName, 'Slc32a1'),i} == 0
    indTbxSlc(i) = 1;
    else
    indTbxSlc(i) = 0;
    end
end

%construct a count matrix for further analysis:
%for just OB neurons -  matr_TbxSlc = count_matr(:,logical(indTbxSlc));
%for all neurons - matr_TbxSlc = count_matr(:,:);

matr_TbxSlc = count_matr(:,:); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------OVERDISPERSED GENES------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find the top of overdispersed genes using the external "selectGenes"
%function  

%the output of this step is saved as "top_overd_Genes_for_all.mat"

top_overd_Genes = selectGenes(matr_TbxSlc,20); %external function
top100overd_Genes = top_overd_Genes(1:100,:);
top500overd_Genes = top_overd_Genes(1:500,:); 
top1000overd_Genes = top_overd_Genes(1:1000,:);


%extract just counts from the tables 
counts = matr_TbxSlc{1:end, 2:end};    
counts_top100 = top100overd_Genes{1:end, 2:(end-7)};
counts_top500 = top500overd_Genes{1:end, 2:(end-7)};
counts_top1000 = top1000overd_Genes{1:end, 2:(end-7)};

%clustegram
HMobj = clustergram(counts_top500','ColumnLabels',top500overd_Genes.(1)',...
        'RowLabels',top500overd_Genes.Properties.VariableNames(2:44),...
        'Symmetric',false,'Colormap',redbluecmap,'DisplayRange',max(max(counts_top500)));
 
sampleInd = 1:(size(counts_top500,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------PCA------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%we use top 500 overdispersed genes for PCA as an input for the
%function "pca"

[pc, zscores, pcvars] = pca(counts_top500');
pcvarsPers = pcvars./sum(pcvars) * 100;
pcvarsPersRound = round(pcvarsPers);
figure();
scatter(zscores(:,1),zscores(:,2),10,'filled','b');
xName = ['PC1 (' num2str(pcvarsPersRound(1)) '%)'];
yName = ['PC2 (' num2str(pcvarsPersRound(2)) '%)'];
%zName = ['PC2 (' num2str(pcvarsPersRound(3)) '%)'];
xlabel(xName);
ylabel(yName);
%zlabel(zName);
title('Principal Component Scatter Plot');
dx = 0.3; dy = 0.3; 
%text(double(zscores(:,1))+dx,double(zscores(:,2))+dy,num2cell(sampleNames));  %sampleInd


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------tSNE------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%we use top 500 overdispersed genes for tSNE as an input for the expertal
%function "tsne"
%Remark: the "tsne" function performs PCA first
%the output of this step is saved as "tsne_final.mat"

tsne_prj = tsne(counts_top500',2,10,33);  %tsne(X, labels, no_dims, initial_dims, perplexity) 
figure();
scatter(tsne_prj(:,1),tsne_prj(:,2),10,'filled','b');
xName = ['tSNE1'];
yName = ['tSNE2'];
xlabel(xName);
ylabel(yName);
title('tSNE:top500');
text(tsne_prj(:,1)+0.1,tsne_prj(:,2)+0.1,num2cell(sampleNames));  %sampleNames sampleInd



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------CLUSTERING------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%perform clustering of the data using the output of the tsne step

clusterer = HDBSCAN(tsne_prj); %create an instance of the HDBSCAN cluster object
% (1) directly set the parameters
clusterer.minpts = 13;           %indicates which nearest neighbor to use for calculating the core-distances of each point in X
clusterer.minclustsize = 5;     %parameter that limits the minimium size of any potential cluster
clusterer.outlierThresh = 0.85; %value between 0 and 1 used to evaluate each point as a potential outlier
clusterer.fit_model(); 			% trains a cluster hierarchy
clusterer.get_best_clusters(); 	% finds the optimal "flat" clustering scheme
clusterer.get_membership();		% assigns cluster labels to the points in X

clusterer.run_hdbscan( clusterer.minpts,clusterer.minclustsize,[],0.85 );  % (2) call run_hdbscan() with optional inputs

figure();
clusterer.plot_tree();
figure();
clusterer.plot_clusters(); 


%%%%%%%%%%%%%%%%%%%%%%%------ANALYSIS OF A CLUSTER------%%%%%%%%%%%%%%%%%%%%%%%

%same analysis on the particular cluster using its index from the
%clustering output

labels = double(clusterer.labels);
clust_numb = 11 %[1, 5, 6, 11];
label_ind = zeros(1,length(labels));
for i = 1:length(labels)
    if  ismember(labels(i), clust_numb);       %labels(i) == clust_numb
        label_ind(i) = 1;
    else label_ind(i) = 0;
    end
end
label_ind = logical(label_ind);

sampleInd_clust_Tbx = sampleInd(label_ind);

[tab_sig_cl_Tbx markers_cl_Tbx]  = signifGenes(counts_top500,top500overd_Genes,...
                                   sampleInd_clust_Tbx,setdiff(sampleInd,sampleInd_clust_Tbx));

%only 1 cluster
sampleInd_ins = sampleInd(label_ind);

%run tsne for the cluster; the output for projection neurons is saved as "tsne_prj_in.mat"
b = counts_top500(:,label_ind);
tsne_prj_in = tsne(b',[],2,10,5);  %tsne(X, labels, no_dims, initial_dims, perplexity) 
figure();
scatter(tsne_prj_in(:,1),tsne_prj_in(:,2),30,'filled','b');
xName = ['tSNE1'];
yName = ['tSNE2'];
xlabel(xName);
ylabel(yName);
title('tSNE:top500');
%text(tsne_prj_in(:,1)+0.1,tsne_prj_in(:,2)+0.1,num2cell(new_int));  %sampleInd_ins
%text(tsne_prj_in(:,1)+0.1,tsne_prj_in(:,2)+0.1,num2cell(sampleInd_ins));  %sampleInd_ins

%plot expression level of some genes in tsne space

expr_lev_plot('A230065H16Rik',top500overd_Genes(:,logical([1,label_ind,1,1,1,1,1,1,1])),tsne_prj_in, []);
expr_lev_plot('Cck',top500overd_Genes(:,logical([1,label_ind,1,1,1,1,1,1,1])),tsne_prj_in);
expr_lev_plot('Tbx21',top500overd_Genes(:,logical([1,label_ind,1,1,1,1,1,1,1])),tsne_prj_in);
expr_lev_plot('Slc17a7',top500overd_Genes(:,logical([1,label_ind,1,1,1,1,1,1,1])),tsne_prj_in);
expr_lev_plot('Slc17a6',top500overd_Genes(:,logical([1,label_ind,1,1,1,1,1,1,1])),tsne_prj_in);
expr_lev_plot('Pkia',top_overd_Genes(:,logical([1,label_ind,1,1,1,1,1,1,1])),tsne_prj_in);
expr_lev_plot('Pkig',top_overd_Genes(:,logical([1,label_ind,1,1,1,1,1,1,1])),tsne_prj_in);

%perform clustering on the tsne output
                               
clusterer = HDBSCAN(tsne_prj_in); %create an instance of the HDBSCAN cluster object
% (1) directly set the parameters
clusterer.minpts = 5;           %indicates which nearest neighbor to use for calculating the core-distances of each point in X
clusterer.minclustsize = 3;     %parameter that limits the minimium size of any potential cluster
clusterer.outlierThresh = 0.85; %value between 0 and 1 used to evaluate each point as a potential outlier
clusterer.fit_model(); 			% trains a cluster hierarchy
clusterer.get_best_clusters(); 	% finds the optimal "flat" clustering scheme
clusterer.get_membership();		% assigns cluster labels to the points in X

clusterer.run_hdbscan( clusterer.minpts,clusterer.minclustsize,[],0.85 );  % (2) call run_hdbscan() with optional inputs

figure();
clusterer.plot_tree();
figure();
clusterer.plot_clusters(); 

%construct the probabilities
figure();
pointsize = 30;
colormap_probab = clusterer.score;
scatter(tsne_prj_in(:,1),tsne_prj_in(:,2), pointsize, colormap_probab,'filled');

                              
labels_ins = double(clusterer.labels);
clust_numb = 2;
label_ind_ins = zeros(1,length(labels_ins));
for i = 1:length(labels_ins)
    if labels_ins(i) == clust_numb
        label_ind_ins(i) = 1;
    else label_ind_ins(i) = 0;
    end
end
label_ind_ins = logical(label_ind_ins);                               
 
cl1 = find(label_ind_ins);%sampleInd_ins(label_ind_ins); 
cl2 = find(label_ind_ins);%sampleInd_ins(label_ind_ins); 
cl3 = find(label_ind_ins);%sampleInd_ins(label_ind_ins);
cl4 = find(label_ind_ins);%sampleInd_ins(label_ind_ins); 
cl5 = find(label_ind_ins);%sampleInd_ins(label_ind_ins); 
cl6 = find(label_ind_ins);%sampleInd_ins(label_ind_ins); 

counts_top500_f = counts_top500(:,label_ind);
label_ind_mod = [1 label_ind 1 1 1 1 1 1 1];
label_ind_mod = logical(label_ind_mod);      
top500overd_Genes_f = top500overd_Genes(:,label_ind_mod);


[tab_sig_cl1 markers_cl1]  = signifGenes(counts_top500_f,top500overd_Genes_f,cl1,[cl2 cl3 cl4 cl5 cl6]);
[tab_sig_cl2 markers_cl2]  = signifGenes(counts_top500_f,top500overd_Genes_f,cl2,[cl1 cl3 cl4 cl5 cl6]);
[tab_sig_cl3 markers_cl3]  = signifGenes(counts_top500_f,top500overd_Genes_f,cl3,[cl1 cl2 cl4 cl5 cl6]);
[tab_sig_cl4 markers_cl4]  = signifGenes(counts_top500_f,top500overd_Genes_f,cl4,[cl1 cl2 cl3 cl5 cl6]);
[tab_sig_cl5 markers_cl5]  = signifGenes(counts_top500_f,top500overd_Genes_f,cl5,[cl1 cl2 cl3 cl4 cl6]);
[tab_sig_cl6 markers_cl6]  = signifGenes(counts_top500_f,top500overd_Genes_f,cl6,[cl1 cl2 cl3 cl4 cl5]);


expr_lev_plot('Pkib',tab_sig_cl4,tsne_prj_in,new_int);
expr_lev_plot('Cck',tab_sig_cl4,tsne_prj_in,new_int);
expr_lev_plot('Cdhr1',top500overd_Genes(:,logical([1,label_ind,1,1,1,1,1,1,1])),tsne_prj_in);
expr_lev_plot('Tbx21',top500overd_Genes(:,logical([1,label_ind,1,1,1,1,1,1,1])),tsne_prj_in);
expr_lev_plot('Slc1a2',top500overd_Genes(:,logical([1,label_ind,1,1,1,1,1,1,1])),tsne_prj_in);
expr_lev_plot('Slc17a7',top500overd_Genes(:,logical([1,label_ind,1,1,1,1,1,1,1])),tsne_prj_in);
expr_lev_plot('Nell2',top500overd_Genes(:,logical([1,label_ind,1,1,1,1,1,1,1])),tsne_prj_in);
expr_lev_plot('Nmb',top500overd_Genes(:,logical([1,label_ind,1,1,1,1,1,1,1])),tsne_prj_in);
expr_lev_plot('Ptn',top500overd_Genes(:,logical([1,label_ind,1,1,1,1,1,1,1])),tsne_prj_in);
