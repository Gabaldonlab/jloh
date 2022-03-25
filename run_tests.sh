#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -c 48

WD="/gpfs/projects/bsc40/current/mschiavi/jloh/results"
cd ${WD}

module load gcc/7.2.0 hisat2/2.1.0 samtools/1.13 bcftools/1.8 htslib/1.8 bedtools/2.29.2


# variables
REF_GENOME="/gpfs/projects/bsc40/current/mschiavi/jloh/raw_data/S_cerevisiae.fa"
unset DIVS; declare -a DIVS=(0.01 0.03 0.05 0.10 0.15 0.20)
unset LOHS; declare -a LOHS=(0.1 0.2 0.3)
unset SUBGENOMES; declare -a SUBGENOMES=(Sc Sd)
unset MINLENS; declare -a MINLENS=(100 250 500 1000)
unset COVS; declare -A COVS=(["10X"]=1000000 ["30X"]=3000000)
unset SNPS; declare -a SNPS=(2 3 4 5)
LOH_MEAN_LEN="5000"

for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          cd ${WD}
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          if [ ! -d ${WD}/${OUT} ]; then mkdir ${WD}/${OUT}; fi
          cd ${WD}/${OUT}
          unset X
          declare -a X=(sim_genome reads variants variants/pileups LOH mapping TP_FP)
          for DIR in ${X[@]}
          do
            if [ ! -d ${DIR} ]; then mkdir ${DIR}; fi
          done
        done
      done
    done
  done
done


# link genome
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          cd ${WD}
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          cd ${WD}/${OUT}/sim_genome
          if [ ! -f Sc.fa ]
          then
            ln -s ${REF_GENOME} Sc.fa
          fi
        done
      done
    done
  done
done


# mutate genome
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          cd ${WD}
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          echo "Mutating genome in ${OUT}..."
          cd ${WD}/${OUT}/sim_genome
          { jloh sim --fasta Sc.fa \
          --threads 40 \
          --out Sd.fa \
          --hybrid \
          --divergence ${DIV} \
          --loh ${LOH} \
          --loh-mean-length ${LOH_MEAN_LEN} \
          --chr-name-mod S_cere,S_dive; } \
          &> Sd.fa.stderr
        done
      done
    done
  done
done

# concatenate the two genomes for IGV viewing
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          cd ${WD}
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          echo "Concatenating the two genomes in ${OUT}..."
          cd ${WD}/${OUT}/sim_genome
          cat Sc.fa Sd.fa > Sc_Sd.fa
          samtools faidx Sc_Sd.fa
        done
      done
    done
  done
done

# simulate reads from both
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          cd ${WD}
          NUM_SIM_READS="${COVS[${COV}]}"
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          echo "Simulating reads in ${OUT}..."
          cd ${WD}/${OUT}/reads
          wgsim -e 0.002 -d 1000 -s 200 -N ${NUM_SIM_READS} -1 125 -2 125 -r 0.0 -R 0.0 \
          -X 0.0 ../sim_genome/Sc_Sd.fa Sc_Sd.R1.fastq Sc_Sd.R2.fastq \
          &> Sc_Sd.wgsim.stderr
        done
      done
    done
  done
done

# index
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          cd ${WD}
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          echo "Indexing genome in ${OUT}..."
          cd ${WD}/${OUT}/sim_genome
          hisat2-build -p 48 Sc.fa Sc.fa &> Sc.fa.stderr &&
          hisat2-build -p 48 Sd.fa Sd.fa &> Sd.fa.stderr &&
          samtools faidx Sc.fa &&
          samtools faidx Sd.fa
        done
      done
    done
  done
done

# map
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          cd ${WD}
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          echo "Mapping reads in ${OUT}..."
          cd ${WD}/${OUT}/mapping
          for SUB in ${SUBGENOMES[@]}
          do
            hisat2-align-s -p 48 --ignore-quals --no-spliced-alignment --score-min L,0.0,-1.0 \
            -I 800 -X 1200 -x ../sim_genome/${SUB}.fa -1 ../reads/Sc_Sd.R1.fastq -2 ../reads/Sc_Sd.R2.fastq \
            -S ${SUB}.sam &> ${SUB}.sam.stderr
          done
        done
      done
    done
  done
done

# filter
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          cd ${WD}
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          echo "Filtering and sorting mapping records in ${OUT}..."
          cd ${WD}/${OUT}/mapping
          for SUB in ${SUBGENOMES[@]}
          do
            samtools view -@ 48 -h -b -F 0x0100 -F 0x4 --output ${SUB}.f.bam ${SUB}.sam &> ${SUB}.samtools.filter.stderr &&
            samtools sort -@ 48 -T ${SUB} -O bam -o ${SUB}.fs.bam ${SUB}.f.bam &> ${SUB}.samtools.sort.stderr
          done
        done
      done
    done
  done
done

# index
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          cd ${WD}
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          echo "Indexing bams for samtools in ${OUT}..."
          cd ${WD}/${OUT}/mapping
          for SUB in ${SUBGENOMES[@]}
          do
            samtools index ${SUB}.fs.bam
          done
        done
      done
    done
  done
done

# pileup
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          cd ${WD}
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          cd ${WD}/${OUT}/variants/pileups
          echo "Performing pileup in ${OUT}..."
          for SUB in ${SUBGENOMES[@]}
          do
            if [ -f ${SUB}.mpileup.vcf.done ]; then rm ${SUB}.mpileup.vcf.done; fi &&
            bcftools mpileup --fasta-ref ../../sim_genome/${SUB}.fa \
            --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
            --output ${SUB}.mpileup.vcf --output-type v --skip-indels \
            ../../mapping/${SUB}.fs.bam &> ${SUB}.mpileup.vcf.stderr &&
            touch ${SUB}.mpileup.vcf.done &
          done
        done
      done
    done
  done
done &&
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          cd ${WD}/${OUT}/variants/pileups
          for SUB in ${SUBGENOMES[@]}
          do
            while [ ! -f ${SUB}.mpileup.vcf.done ]
            do
              sleep 1
            done
          done
        done
      done
    done
  done
done

# call
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          cd ${WD}
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          cd ${WD}/${OUT}/variants
          echo "Calling variants in ${OUT}..."
          for SUB in ${SUBGENOMES[@]}
          do
            bcftools call --threads 4 --multiallelic-caller --variants-only \
            --output ${SUB}.raw.vcf --output-type v pileups/${SUB}.mpileup.vcf \
            &> ${SUB}.call.raw.vcf.stderr &&
            touch ${SUB}.call.raw.vcf.done &
          done
        done
      done
    done
  done
done &&
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          cd ${WD}/${OUT}/variants
          for SUB in ${SUBGENOMES[@]}; do
            while [ ! -f ${SUB}.call.raw.vcf.done ]
            do
              sleep 1
            done
          done
        done
      done
    done
  done
done

# filter
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          cd ${WD}
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          cd ${WD}/${OUT}/variants
          echo "Filtering variants in ${OUT}..."
          for SUB in ${SUBGENOMES[@]}
          do
            /gpfs/projects/bsc40/mschiavi/software/all2vcf/all2vcf filter_vcf \
            --input-file ${SUB}.raw.vcf --map-qual-zero-frac 0.05 \
            --output-file ${SUB}.f.vcf --quality 20 --alt-frac 0.05 --min-depth 4 \
            &> ${SUB}.all2vcf_filter.stderr &&
            /gpfs/projects/bsc40/mschiavi/software/all2vcf/all2vcf frequency \
            --in ${SUB}.f.vcf --out ${SUB}.ff.vcf \
            &> ${SUB}.all2vcf_filter.freq.stderr
          done
        done
      done
    done
  done
done

# JLOH g2g
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          cd ${WD}
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          cd ${WD}/${OUT}/LOH
          echo "Aligning the two genomes to find SNPs in ${OUT}..."
          /gpfs/projects/bsc40/mschiavi/software/jloh/source/jloh g2g \
          --sensitive \
          --ref-A ../sim_genome/Sc.fa --ref-B ../sim_genome/Sd.fa --est-divergence ${DIV} \
          --out . \
          --min-length ${MINLEN} \
          &> Sc_Sd.regions.bed.stderr &&
          cat A.bed B.bed > Sc_Sd.regions.bed &&
          touch Sc_Sd.regions.bed.done &
        done
      done
    done
  done
done &&
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          cd ${WD}/${OUT}/LOH
          while [ ! -f Sc_Sd.regions.bed.done ]
          do
            sleep 1
          done
        done
      done
    done
  done
done

# JLOH extract
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          cd ${WD}
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          cd ${WD}/${OUT}/LOH
          echo "Extracting LOH blocks in ${OUT}..."
          /gpfs/projects/bsc40/mschiavi/software/jloh/source/jloh extract \
          --hybrid --threads 48 --vcfs ../variants/Sc.ff.vcf ../variants/Sd.ff.vcf \
          --bams ../mapping/Sc.fs.bam ../mapping/Sd.fs.bam \
          --refs ../sim_genome/Sc.fa ../sim_genome/Sd.fa \
          --sample Sc_Sd --output-dir . --regions Sc_Sd.regions.bed \
          --min-length ${MINLEN} --min-snps 2 \
          &> LOH_extraction.stderr &&
          { cat Sc_Sd.LOH_blocks.A.tsv; tail -n+2 Sc_Sd.LOH_blocks.B.tsv; } \
          > Sc_Sd.LOH_blocks.tsv &&
          cat Sc_Sd.LOH_blocks.A.bed Sc_Sd.LOH_blocks.B.bed \
          > Sc_Sd.LOH_blocks.bed &&
          cat Sc_Sd.LOH_candidates.A.bed Sc_Sd.LOH_candidates.B.bed > Sc_Sd.LOH_candidates.bed
        done
      done
    done
  done
done

# remove tmp folders and tmp files that are useless at this point
echo "Removing tmp files and folders..."
cd $WD
find . -type d -name "tmp_*" | xargs rm -r
find . -type d -name "pileups" | xargs rm -r
find . -type f -name "*sam" | xargs rm
find . -type f -name "*.f.bam" | xargs rm
find . -type d -name "reads" | xargs rm -r
find . -type f -name "*raw.vcf" | xargs rm
find . -type f -name "*.f.vcf" | xargs rm

# TP / FP verification
for DIV in ${DIVS[@]}
do
  for LOH in ${LOHS[@]}
  do
    for MINLEN in ${MINLENS[@]}
    do
      for COV in ${!COVS[@]}
      do
        for MIN_SNPS in ${SNPS[@]}
        do
          cd ${WD}
          OUT="snps_${MIN_SNPS}.minlen_${MINLEN}.cov_${COV}.run.div_${DIV}.loh_${LOH}"
          cd ${WD}/${OUT}/TP_FP
          echo "Calculating precision and recall in ${OUT}..."
          cat ../LOH/Sc_Sd.LOH_blocks.bed | awk '$1 ~ /S_dive/' > selected.bed &&
          cat ../LOH/Sc_Sd.LOH_candidates.bed | awk '$1 ~ /S_dive/' > predicted.bed &&
          { cat ../sim_genome/Sd.fa.non_divergent | awk -v x=${MINLEN} '$3-$2 >= x'; \
          cat ../sim_genome/Sd.fa.lohs | awk -v x=${MINLEN} '$3-$2 >= x'; } | \
          sort -k1,1 -k2n,2 -k3nr,3 > generated.bed &&
          python3 /gpfs/projects/bsc40/current/mschiavi/jloh/scripts/src/tp_fp_rates.py \
          --generated generated.bed \
          --predicted predicted.bed \
          --selected selected.bed \
          --out stats.tsv \
          &> calculations.stderr
        done
      done
    done
  done
done

# make summary
cd ${WD} &&
{ \
echo -e "MinSnps\tMinLength\tCoverage\tDivergence\tLOH\tTP\tFP\tFN\tPrecision\tRecall" &&
for i in $(find . -type f -name "stats.tsv")
do
  echo -e "${i}\t$(cat $i | grep "predicted")"
done | \
sed -e 's/\/TP_FP\/stats.tsv//' | \
sed -e 's/\.\///' | \
sed -e 's/snps_//' | \
sed -e 's/.minlen_//' | \
sed -e 's/.run.div_/\t/' | \
sed -e 's/.loh_/\t/' | \
sed -e 's/.cov_/\t/' | \
sed -e 's/predicted//' | \
sort -k2g,2 -k3g,3; } \
> TP_FP_stats.tsv
