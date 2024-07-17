function [subsets_unique, pfault_inst_unique, pfault_exp_unique, idx_a] = find_unique_subsets(subsets, pfault_inst, pfault_exp)

[subsets_unique, idx_a, idx_c] = unique(subsets,'rows','stable');
pfault_inst_unique = zeros(length(idx_a),1);
pfault_exp_unique = zeros(length(idx_a),1);


%subsets = subsets_unique(idx_c)
for i =1:max(idx_c)
    id=(idx_c==i);
    pfault_inst_unique(i) = sum(pfault_inst(id));
    pfault_exp_unique(i) = sum(pfault_exp(id));   
end
