function top_overd_Genes = selectGenes(matr_TbxSlc,numb_bins)

mtr = matr_TbxSlc{:,2:end};

sz_matr = size(mtr);
%remove strings with low expression level
for i = 1:sz_matr(1)
    if sum(mtr(i,2:end)) == 0  %< 0.1      %1
        ind_del(i) = 0;
    else 
         ind_del(i) = 1;
    end
    disp({'Iteration (removing)' num2str(i) ' of ' num2str(sz_matr(1))});
end

matr_TbxSlc_nonZer = mtr(logical(ind_del'),:);
tab_TbxSlc_nonZer = matr_TbxSlc(logical(ind_del'),:);

sz = size(matr_TbxSlc_nonZer);
for i = 1: sz(1)
meanExpr(i) = mean(matr_TbxSlc_nonZer(i,:));
varExpr(i) = var(matr_TbxSlc_nonZer(i,:));
FanoFactor(i) = log10(varExpr(i) / meanExpr(i));
disp({'Iteration (statist. values)' num2str(i) ' of ' num2str(sz_matr(1))});
end

tab_TbxSlc_nonZer.meanExpr = meanExpr';
tab_TbxSlc_nonZer.varExpr = varExpr';
tab_TbxSlc_nonZer.FanoFactor = FanoFactor';


step = (max(meanExpr) - min(meanExpr))/numb_bins;
for i = 1:(numb_bins+1)
   bins(i) = min(meanExpr) + step*(i-1);
   disp({'Iteration (bins)' num2str(i) ' of ' num2str(sz_matr(1))});
end

for i = 1:numb_bins
     for j = 1:sz(1)
         if meanExpr(j) >= bins(i) && meanExpr(j) <= bins(i+1)
             %matr_TbxSlc_nonZer.binInd(j) = i;  
             binInd(j) = i; 
         end
     end
     disp({'Iteration (bin_ind)' num2str(i) ' of ' num2str(sz_matr(1))});
end

tab_TbxSlc_nonZer.binInd = binInd';

for i = 1:numb_bins
     for j = 1:sz(1)
         if binInd(j) ==i
             Ind(j,i) = 1;
         else Ind(j,i) = 0;
     end
     end
     disp({'Iteration (Ind)' num2str(i) ' of ' num2str(sz_matr(1))});
end

for i = 1:numb_bins
  indBin = logical(Ind(:,i));
  mean_FrfactorBin(i) = mean(FanoFactor(indBin)');
  std_FrfactorBin(i) = std(FanoFactor(indBin)');
end

for j = 1:sz(1)
    for i = 1:numb_bins
       if  binInd(j) == i
           mean_FrfactorBin(j) = mean_FrfactorBin(i);
           std_FrfactorBin(j) = std_FrfactorBin(i);
       end
    end
end 

tab_TbxSlc_nonZer.binInd = binInd';
tab_TbxSlc_nonZer.meanFrfactorBin = mean_FrfactorBin';
tab_TbxSlc_nonZer.std_FrfactorBin = std_FrfactorBin';

for i = 1: sz(1)
   d(i) = (FanoFactor(i) -  mean_FrfactorBin(i))/ std_FrfactorBin(i);
  disp({'Iteration (final)' num2str(i) ' of ' num2str(sz_matr(1))});
end

tab_TbxSlc_nonZer.d = d';

disp('Final step...');

[out,idx] = sort(d, 'descend');

idx_f = idx(1:1000);

top_overd_Genes = tab_TbxSlc_nonZer(idx_f,:);












%{

sz_matr = size(matr_TbxSlc);
%remove strings with low expression level
for i = 1:sz_matr(1)
    if mean(matr_TbxSlc{i,2:end}) == 0  %< 0.1      %1
        ind_del(i) = 0;
    else 
         ind_del(i) = 1;
    end
    disp({'Iteration (removing)' num2str(i) ' of ' num2str(sz_matr(1))});
end

matr_TbxSlc_nonZer = matr_TbxSlc(logical(ind_del'),:);


%remove zero strings
%counts_prjNeur = matr_TbxSlc{1:end,2:end};
%matr_TbxSlc_nonZer = matr_TbxSlc;
%matr_TbxSlc_nonZer( ~any(counts_prjNeur,2), : ) = [];
%
%TPM_resc_prj = log2(counts_prj+1);

%counts_TbxSlc_nonZer = matr_TbxSlc_nonZer{:,2:end};

sz = size(matr_TbxSlc_nonZer);
for i = 1: sz(1)
meanExpr(i) = mean(matr_TbxSlc_nonZer{i,2:end});
varExpr(i) = var(matr_TbxSlc_nonZer{i,2:end});
FanoFactor(i) = log10(varExpr(i) / meanExpr(i));
disp({'Iteration (statist. values)' num2str(i) ' of ' num2str(sz_matr(1))});
end

matr_TbxSlc_nonZer.meanExpr = meanExpr';
matr_TbxSlc_nonZer.varExpr = varExpr';
matr_TbxSlc_nonZer.FanoFactor = FanoFactor';


step = (max(meanExpr) - min(meanExpr))/numb_bins;
for i = 1:(numb_bins+1)
   bins(i) = min(meanExpr) + step*(i-1);
   disp({'Iteration (bins)' num2str(i) ' of ' num2str(sz_matr(1))});
end

for i = 1:numb_bins
     for j = 1:sz(1)
         if matr_TbxSlc_nonZer.meanExpr(j) >= bins(i) && matr_TbxSlc_nonZer.meanExpr(j) <= bins(i+1)
             %matr_TbxSlc_nonZer.binInd(j) = i;  
             binInd(j) = i; 
         end
     end
     disp({'Iteration (bin_ind)' num2str(i) ' of ' num2str(sz_matr(1))});
end

matr_TbxSlc_nonZer.binInd = binInd';

for i = 1:numb_bins
     for j = 1:sz(1)
         if matr_TbxSlc_nonZer.binInd(j) ==i
             Ind(j,i) = 1;
         else Ind(j,i) = 0;
     end
     end
     disp({'Iteration (Ind)' num2str(i) ' of ' num2str(sz_matr(1))});
end

for i = 1:numb_bins
  indBin = logical(Ind(:,i));
  mean_FrfactorBin(i) = mean(matr_TbxSlc_nonZer.FanoFactor(indBin)',2);
  std_FrfactorBin(i) = std(matr_TbxSlc_nonZer.FanoFactor(indBin)');
end


for j = 1:sz(1)
    for i = 1:numb_bins
       if  matr_TbxSlc_nonZer.binInd(j) == i
           matr_TbxSlc_nonZer.meanFrfactorBin(j) = mean_FrfactorBin(i);
           matr_TbxSlc_nonZer.stdFrfactorBin(j) = std_FrfactorBin(i);
       end
    end
end 

for i = 1: sz(1)
 matr_TbxSlc_nonZer.d(i) = (matr_TbxSlc_nonZer.FanoFactor(i) -  matr_TbxSlc_nonZer.meanFrfactorBin(i))/ matr_TbxSlc_nonZer.stdFrfactorBin(i);
  disp({'Iteration (final)' num2str(i) ' of ' num2str(sz_matr(1))});
end

top_overd_Genes = sortrows(matr_TbxSlc_nonZer,'d','descend');

%}

%top500overd_Genes_ID = top500overd_Genes.Properties.RowNames;


end