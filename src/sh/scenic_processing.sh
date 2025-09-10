### pyscenic pipeline
## inputs
input_loom=result/Regulon/sample.loom

## outputs
grn_output=result/Regulon/out_adj.tsv
ctx_output=result/Regulon/out_reg.tsv
loom_output=result/Regulon/out_SCENIC.loom

## reference
f_tfs=../cisTarget_db/mus_mm10_tfs.motifs-v10.txt
f_motif_path=../cisTarget_db/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
f_db_500bp=../cisTarget_db/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
f_db_10kb=../cisTarget_db/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather


#1.1 grn
pyscenic grn \
 --num_workers 30 \
 --output $grn_output \
 --method grnboost2 \
 $input_loom \
 $f_tfs

#1.2 cistarget
pyscenic ctx \
 $grn_output \
 $f_db_500bp $f_db_10kb \
 --annotations_fname $f_motif_path \
 --expression_mtx_fname $input_loom \
 --output $ctx_output \
 --num_workers 30

#1.3 AUCell
pyscenic aucell \
 $input_loom \
 $ctx_output \
 --output $loom_output \
 --num_workers 30