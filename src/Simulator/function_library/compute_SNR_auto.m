function cut_off=compute_SNR_auto(S_matrix)

log10_s_matrix=log10(S_matrix);

max_cutoff=max(log10_s_matrix(end,:));
min_cutoff=min(min(log10_s_matrix(3,:)),max(log10_s_matrix(:))-2);

cutoff_candidates=linspace(min_cutoff,max_cutoff,100);
for cand=1:length(cutoff_candidates)
    for s_no=1:size(S_matrix,2)
        s_noise(cand,s_no)=find(log10_s_matrix(:,s_no)<=cutoff_candidates(cand),1,'first')+sum(isinf(log10_s_matrix(:,s_no)));
%         sum(isinf(log10_s_matrix(:,s_no)))
    end
    span(cand)=max(s_noise(cand,:))-min(s_noise(cand,:));
end
[~,cand_cut]=min(span);
cut_off=cutoff_candidates(cand_cut);
