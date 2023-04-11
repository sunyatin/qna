MODEL=$1
ITER=$2
OVERWRITE=$3
# either 0 or 1
MAX_CHROM=20
TEMPDIR=tmpdir/${LOGNAME}/20x30Mbp.1_2e-8
#${SLURM_JOBID}
OUTDIR=/users/p23002/p23002tr/qian/Final.Blake/data/20x30Mbp.1_2e-8/${MODEL}
#___
echo "===================================="
echo $MODEL
echo $ITER
echo $OVERWRITE
#___

mkdir -p $TEMPDIR/$MODEL
mkdir -p $OUTDIR
#___
CORE="$SCRIPT_SIMULATE \
-o $TEMPDIR/$MODEL/${MODEL}_${ITER} \
--write_geno1 \
--bin_dir $BIN_DIR \
-i @published/1.2e-8.est \
--mut_model BinaryMutationModel(False) \
--algo hudson \
--seed $ITER \
--model $MODEL \
-g 20 30"
##
if [ "$MODEL" == "kamm19" ]; then
SAMPLES="Mbuti:50:0:YRI Sardinian:50:0:CEU Neanderthal:1:50000:Vindija"
elif [ "$MODEL" == "ragsdale19" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Neanderthal:1:50000:Vindija"
elif [ "$MODEL" == "jacobs19" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU NeaA:1:50000:Vindija"
elif [ "$MODEL" == "gower21" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Nea:1:50000:Vindija"
elif [ "$MODEL" == "iasi21" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU NEA:1:50:Vindija"
elif [ "$MODEL" == "durvasula20" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Neanderthal:1:50000:Vindija"
elif [ "$MODEL" == "fu14" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Nea_1:1:50000:Vindija"
elif [ "$MODEL" == "yang12" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU NEA:1:50000:Vindija"
elif [ "$MODEL" == "skov20" ]; then
SAMPLES="Africa:50:0:YRI Iceland:50:0:CEU Vindija:1:50000:Vindija"
elif [ "$MODEL" == "moorjani16" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Nea_1:1:50000:Vindija"
elif [ "$MODEL" == "schaefer21" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Vindija:1:50000:Vindija"
else
echo "$MODEL not implemented"
fi
##
CMD="$CORE --samples $SAMPLES"
##
if [ "$OVERWRITE" -eq "0" ]; then
	if [ -f $TEMPDIR/$MODEL/${MODEL}_${ITER}.snp.gz ]; then
		LAST_CHR=$(zcat $TEMPDIR/$MODEL/${MODEL}_${ITER}.snp.gz | tail -n1 | cut -f2 -d' ')
		if [ $LAST_CHR -eq $MAX_CHROM ]; then
			echo "PASS"
		else
			echo "SIMULATE"
			python3.7 $CMD
		fi
	else
		echo "SIMULATE"
		python3.7 $CMD
	fi
elif [ "$OVERWRITE" -eq "1" ]; then
	echo "SIMULATE"
	python3.7 $CMD
else
	echo "ERROR IN 3rd ARGUMENT"
fi
#___
cd $TEMPDIR/$MODEL
CMD="$SCRIPT_SUMSTATS \
--bin_dir $BIN_DIR \
-p /home/rtournebize/gitdir/qna/param_files/stats_MAC_diplo.spar \
--stats --SE --sprime --crf --psmc --ld \
--ceu CEU --yri YRI --vindija Vindija \
--psmc_pops CEU YRI \
-i ${MODEL}_${ITER}"
##
if [ "$OVERWRITE" -eq "0" ]; then
	if [ -f $TEMPDIR/$MODEL/${MODEL}_${ITER}.crf.stats ]; then
		echo "PASS"
	else
		echo "SUMSTATS"
		python3.7 $CMD
		cp ${MODEL}_${ITER}.*stats ${MODEL}_${ITER}.par $OUTDIR
	fi
elif [ "$OVERWRITE" -eq "1" ]; then
	echo "SUMSTATS"
	python3.7 $CMD
	cp ${MODEL}_${ITER}.*stats ${MODEL}_${ITER}.par $OUTDIR
else
	echo "ERROR IN 3rd ARGUMENT"
fi
#___
