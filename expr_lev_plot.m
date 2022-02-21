function y = expr_lev_plot(gene,tab_sig,tsne_prj,sampleNames)

gene_expr_level_tab = tab_sig{ismember(tab_sig.GeneName, gene),2:44};
figure();
pointsize = 70;
colormap_gene = gene_expr_level_tab;
scatter(tsne_prj(:,1),tsne_prj(:,2), pointsize, colormap_gene,'filled');
colorbar
title([gene ' expression level']);
xName = ['tSNE1'];
yName = ['tSNE2'];
xlabel(xName);
ylabel(yName);
text(tsne_prj(:,1),tsne_prj(:,2),sampleNames);

end