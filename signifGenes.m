function [tab_sig markers] = signifGenes(counts_top500,top500overd_Genes,ind_cl1_p,ind_cl2_p)

for i = 1:size(counts_top500,1)
pVal(i) = ranksum(counts_top500(i,ind_cl1_p),counts_top500(i,ind_cl2_p));
disp({'Iteration (pVal)' num2str(i) ' of ' size(counts_top500,1)});
end
padj = mafdr(pVal,'BHFDR',true);
top500overd_Genes.padj = padj';

sig = top500overd_Genes.padj  < 0.05 ;

tab_sig = top500overd_Genes;%top500overd_Genes(sig,:);

%calculate median expression within cluster of interest
for i = 1:size(tab_sig,1)
    med_expr_cl1(i) = median(tab_sig{i,ind_cl1_p+1});
    med_expr_cl2(i) = median(tab_sig{i,ind_cl2_p+1});
    disp({'Iteration (med expr)' num2str(i) ' of ' size(tab_sig,1)});
end

tab_sig.medExpr_1 =  med_expr_cl1';
tab_sig.medExpr_2 =  med_expr_cl2';

threshold = 3;
%filter for clusters
for i = 1:size(tab_sig,1)
    if tab_sig.medExpr_1(i) > threshold
        tab_sig.signEnr_cl1(i) = 1;
    else tab_sig.signEnr_cl1(i) = 0;
    end
    
    if tab_sig.medExpr_2(i) > threshold
        tab_sig.signEnr_cl2(i) = 1;
    else tab_sig.signEnr_cl2(i) = 0;
    end
    disp({'Iteration (filter)' num2str(i) ' of ' size(tab_sig,1)});
end

%fraction of cells expr the gene for each cluster
for i = 1:size(tab_sig,1)
    
      tab_sig.FractCells_cl1(i) =  sum(tab_sig{i,ind_cl1_p+1} > threshold) / length(ind_cl1_p);
      tab_sig.FractCells_cl2(i) = sum(tab_sig{i,ind_cl2_p+1} > threshold) / length(ind_cl2_p);
      
    if  tab_sig.FractCells_cl1(i) >= 0.5 && tab_sig.FractCells_cl2(i) <= 0.1
        tab_sig.Marker_cl1(i) = 1;
    else tab_sig.Marker_cl1(i) = 0;
    end
    
    if  tab_sig.FractCells_cl2(i) >= 0.5 && tab_sig.FractCells_cl1(i) <= 0.1
        tab_sig.Marker_cl2(i) = 1;
    else tab_sig.Marker_cl2(i) = 0;
    end
    disp({'Iteration (final)' num2str(i) ' of ' size(tab_sig,1)});
end

markers = tab_sig.GeneName(logical(tab_sig.Marker_cl1));

end