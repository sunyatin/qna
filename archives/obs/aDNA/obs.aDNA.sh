# To run the following script, the AADR v54.2 1240K dataset must have been downloaded and
# converted to EIGENSTRAT (non binarized/compressed *.geno file) format using convertf.

# Then, overwrite the `v54.1_1240K_public.ind` file with the `v54.1_1240K_public_MODIFIED.ind` file
# from this directory (keeping the name `v54.1_1240K_public.ind`)
# This *.ind file has a few population labels changed to avoid analyzing duplicated individual entries.

# Then, change the PREFIX_AADR_DATASET path accordingly:
PREFIX_AADR_DATASET=/home/rtournebize/data/v54.1.1240K/unbinarized/v54.1_1240K_public

for IND in Russia_Yana_old2_UP.SG Russia_Kostenki14 Romania_Oase_UP_enhanced Russia_Ust_Ishim.DG French.DG Luxembourg_Loschbour.DG Germany_EN_LBK_Stuttgart.DG Russia_MA1_HG.SG Bulgaria_BachoKiro_LatePleistocene_BB7 Bulgaria_BachoKiro_LatePleistocene_CC7335
do
python3.7 scripts/sumstats.py \
--bin_dir bin \
-p param_files/further_stats.spar \
--stats --SE --f4 --ld1 \
--ceu $IND --yri Yoruba.DG --vindija Vindija_Neanderthal.DG --altai Altai_Neanderthal.DG --ancestral Ancestor.REF \
-i $PREFIX_AADR_DATASET \
-o archives/obs/stats_aDNA/$IND \
--chrom 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 \
--na_rm | tee archives/obs/stats_aDNA/$IND.log
done

#___
