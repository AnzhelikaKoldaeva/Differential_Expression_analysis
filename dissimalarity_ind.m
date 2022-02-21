bin_matrix = readtable( 'ECFP_1024_m0-2_b1_c.txt' ) ;

for i = 1:size(bin_matrix,1)
    for j = 1:size(bin_matrix,1)
        and_comp = and(bin_matrix{i,2:end}, bin_matrix{j,2:end});
        or_comp = or(bin_matrix{i,2:end}, bin_matrix{j,2:end});
        similar_ind_matrix(i,j) = sum(and_comp) / sum(or_comp);
    end
end

var_names = {'Acetophen', 'EthylButyr', 'EthylPropionate', 'Eugenol',...
             'IsopropAcet', 'MethylAnthran','MethylBenz','MethylButyr',...
             'MethylSalicyl','MethylTigl','Salisylal', 'Vanilin'};
         
 bin_matrix.Properties.RowNames =   var_names;
 
 similar_ind_tab = array2table(similar_ind_matrix);
 similar_ind_tab.Properties.RowNames =   var_names;
 similar_ind_tab.Properties.VariableNames =   var_names;
 
 
 HMobj = clustergram(similar_ind_matrix','ColumnLabels',var_names,...
        'RowLabels',var_names,...
        'Symmetric',false,'Colormap',winter,'DisplayRange',max(max(similar_ind_matrix)));
 