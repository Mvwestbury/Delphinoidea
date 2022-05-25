# Delphinoidea
Commands, pipeline, and scripts for A genomic assessment of the marine-speciation paradox within the toothed whale superfamily Delphinoidea

For Github

## Sequencing read filtering and mapping
- Included script: Modern_mapping_mem_PE.sh

## Finding sex-linked scaffolds/autosomes
- Use satsuma synteny
`satsuma - SatsumaSynteny -m 1 -n 30 -q Referencegenome.fasta -t X_Y_chromosomes.fasta -o output_directory`
- Obtain scaffold names
`cut -f 4 satsuma_summary.chained.out | sort | uniq | sed 's/1_/1\t/g' | cut -f 1 > Reference_XY.txt`
- Create text file with only autosomes and scaffolds >100kb
`samtools faidx Referencegenome.fasta`
`grep -v -f Reference_XY.txt Referencegenome.fasta.fai | awk '$2 >100000 {print $1}' > Reference_noXY_100kb.txt`

## Fasta call
`angsd -dofasta 2 -mininddepth 5 -minmapq 30 -minq 30 -uniqueonly 1 -only_proper_pairs 1 -docounts 1 -rf Reference_noXY_100kb.txt -i bamfile.bam -out bamfile_autosomes_100kb`

## Sliding window trees
- Filtering and RAxML
- Included script: windowTrees.pl but this script has been modified for easier use and is available at https://github.com/achimklittich/WindowTrees
- Filtering based on GC content: countGC.pl
- QuiBL
- ASTRAL
`java -jar astral.5.7.8.jar -i out1k.trs -o out1k_astral.tre`

## Dstatistics-ANGSD
`angsd -doabbababa 2 -nthreads 10 -mininddepth 5 -minmapq 30 -minq 30 -uniqueonly 1 -only_proper_pairs 1 -docounts 1 -b bamlist.txt -rf Lipotes_vexillifer_noXY_1Mb.txt -ref GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna -anc GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna -out Toothed_whales_Dstats -blocksize 1000000`
`Rscript jackKnife.R file=Toothed_whales_Dstats.abbababa indNames=names.txt outfile=Toothed_whales_Dstats.jackknife`

## Dfoil
- Prepare fasta files for Dfoil
- Round down to the nearest 100kb
`cat Lipotes_vexillifer_noXY_100kb.fasta.fai | awk '{x=int($2+0.5);sub(".....$","00000",x);print $1"\t0\t"x}'`

- Extract scaffolds while removing ends of scaffolds
`bedtools getfasta -fi Lipotes_vexillifer_noXY_100kb.fasta -bed Lipotes_vexillifer_noXY_100kb.rounddown.bed > Lipotes_vexillifer_noXY_100kb.rounddown.fasta`

Remove scaffold headers and add new header for the sequence
`sed '/^>/ d' Lipotes_vexillifer_noXY_100kb.rounddown.fasta | sed -e '1i\>Baiji' > labeled.fa`

`bedtools getfasta -fi $file -bed Lipotes_vexillifer_noXY_100kb.rounddown.bed > ${bn}.rounddown.fasta`
`python mvftools.py ConvertFasta2MVF --fasta Fasta_files/All_relabeled_combined.fasta --out All_relabeled_combined.mvf`

`python /groups/hologenomics/westbury/Software/mvftools/mvftools.py CalcPatternCount --mvf /groups/hologenomics/westbury/data/Toothed_whale_Geneflow/Analyses/Dfoil/All_relabeled_combined.mvf --out Finless_Harbour_porpoise_Pilotwhale_Orca_Baiji.100kb.input --windowsize 100000 --sample-labels Finless,Harbour_porpoise,Pilotwhale,Orca,Baiji`

`python /groups/hologenomics/westbury/Software/mvftools/mvftools.py CalcPatternCount --mvf /groups/hologenomics/westbury/data/Toothed_whale_Geneflow/Analyses/Dfoil/All_relabeled_combined.mvf --out Finless_Harbour_porpoise_Beluga_Narwhal.100kb.input --windowsize 100000 --sample-labels Finless,Harbour_porpoise,Beluga,Narwhal,Baiji`
`dfoil.py --infile $file --out ${bn}.100kb.output`
`dfoil_analyze.py $file > $file.result`

## Dsuite
`reformat.sh minlength=100000 in=GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna out=GCF_000442215.1_Lipotes_vexillifer_v1_genomic_100kb.fa`
`wgsim -1 150 -2 150 -r 0 -R 0 -X 0 -e 0 -N 100000000 GCF_000442215.1_Lipotes_vexillifer_v1_genomic_100kb.fa reference_1.fastq reference_2.fastq`
- Make vcf
`bcftools mpileup --threads 5 -f ../../Reference/GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna -Ou -R Lipotes_100kb.txt -b bamlist.txt | bcftools call -mv -Ov --threads 5 -o All.vcf`
`Dsuite Dtrios -t Tree.txt -o testDtri All_test.vcf Sets.txt`
`Dsuite Fbranch Tree.txt testDtri_tree.txt > fbranch.txt`
`dtools.py fbranch.txt Tree.txt`


## Mutation rate
`angsd -doIBS 2 -nthreads 5 -mininddepth 5 -minmapq 30 -minq 30 -uniqueonly 1 -only_proper_pairs 1 -docounts 1 -minInd 9 -makematrix 1 -b ../Dstats/bamlist.txt -rf /groups/hologenomics/westbury/data/Toothed_whale_Geneflow/Reference/Lipotes_vexillifer_noXY_100kb.txt -ref /groups/hologenomics/westbury/data/Toothed_whale_Geneflow/Reference/GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna -out Toothed_whales_IBS` 

## hPSMC
- Merge 2x consensuses into psmcfa file
`/groups/hologenomics/westbury/Software/hPSMC/psmcfa_from_2_fastas.py -b10 -m5 $1_autosomes_100kb.fa $2_autosomes_100kb.fa > $1_$2_hPSMC.psmcfa`
- Run psmc
`psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o $1_$2_hPSMC.psmc $1_$2_hPSMC.psmcfa`
- Get predivergence Ne from the plot
`psmc_plot.pl -g 20 -u 7.90E-09 -R -s 10 test Beluga_Finless_hPSMC.psmc`
- run simulations
`hPSMC_quantify_split_time.py -N 20000 -l 500000 -u 5000000 -p 11 -s 11 -o output_prefix`

## Relative dating
`Use scripts found in whaleDatingScripts directory

## Heterozygosity
`angsd -GL 1 -mininddepth 5 -minmapq 30 -minq 30 -uniqueonly 1 -only_proper_pairs 1 -docounts 1 -i bamfile -ref GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna -P 5 -out $line -doSaf 1 -anc GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna -rf Lipotes_vexillifer_noXY_100kb.txt -baq 1
`realSFS -nSites 20000000 $line.saf.idx -P 10 -tole 1e-8 2> log > $line_20Mb.sfs`

## PSMC
`samtools mpileup -Q 30 -uf reference.fasta bamfile.bam | bcftools call -c - | vcfutils.pl vcf2fq -d 10 | gzip > file_diploid.fq.gz`
`utils/splitfa diploid.psmcfa > diploid.split.psmcfa`
# Create PSMC data with whole sequence data
`psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o diploid.psmc diploid.psmcfa``
# Bootstrap replicates (with multiple threads)
`seq 100 | xargs -P 10 -i ./psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o diploid.round-{}.psmc diploid.split.psmcfa
	# -P = number of threads
	# seq = number of bootstraps
	# -o = outputfile
#Combine the psmc output files
`cat diploid.psmc diploid.round-*.psmc > combined.diploid.psmc`




