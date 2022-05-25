# Delphinoidea
Commands, pipeline, and scripts for A genomic assessment of the marine-speciation paradox within the toothed whale superfamily Delphinoidea.

Note: This is not meant to be standalone commands but examples to be used for recreating the analyses in the manuscript

## Sequencing read filtering and mapping
- Included script: Modern_mapping_mem_PE.sh

## Finding sex-linked scaffolds/autosomes https://github.com/bioinfologics/satsuma2
- Use satsuma synteny

`satsuma - SatsumaSynteny -m 1 -n 30 -q Referencegenome.fasta -t X_Y_chromosomes.fasta -o output_directory`
- Obtain scaffold names

`cut -f 4 satsuma_summary.chained.out | sort | uniq | sed 's/1_/1\t/g' | cut -f 1 > Reference_XY.txt`
- Create text file with only autosomes and scaffolds >100kb

`samtools faidx Referencegenome.fasta`

`grep -v -f Reference_XY.txt Referencegenome.fasta.fai | awk '$2 >100000 {print $1}' > Reference_noXY_100kb.txt`

## Fasta call https://github.com/ANGSD/angsd

`angsd -dofasta 2 -mininddepth 5 -minmapq 30 -minq 30 -uniqueonly 1 -only_proper_pairs 1 -docounts 1 -rf Reference_noXY_100kb.txt -i bamfile.bam -out bamfile_autosomes_100kb`

## Sliding window trees
- Filtering and RAxML
- Included script: windowTrees.pl but this script has been modified for easier use and is available at https://github.com/achimklittich/WindowTrees
- Filtering based on GC content: countGC.pl
- QuiBL (config included)

`QuIBL.py myInputFile2.txt`
- ASTRAL

`java -jar astral.5.7.8.jar -i out1k.trs -o out1k_astral.tre`

## Dstatistics https://github.com/ANGSD/angsd

`angsd -doabbababa 2 -nthreads 10 -mininddepth 5 -minmapq 30 -minq 30 -uniqueonly 1 -only_proper_pairs 1 -docounts 1 -b bamlist.txt -rf Lipotes_vexillifer_noXY_1Mb.txt -ref GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna -anc GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna -out Toothed_whales_Dstats -blocksize 1000000`

`Rscript jackKnife.R file=Toothed_whales_Dstats.abbababa indNames=names.txt outfile=Toothed_whales_Dstats.jackknife`

## Dfoil https://github.com/jbpease/dfoil 
- Prepare fasta files for Dfoil
- Round down to the nearest 100kb

`cat Lipotes_vexillifer_noXY_100kb.fasta.fai | awk '{x=int($2+0.5);sub(".....$","00000",x);print $1"\t0\t"x}'`
- Extract scaffolds while removing ends of scaffolds

`bedtools getfasta -fi Lipotes_vexillifer_noXY_100kb.fasta -bed Lipotes_vexillifer_noXY_100kb.rounddown.bed > Lipotes_vexillifer_noXY_100kb.rounddown.fasta`
Remove scaffold headers and add new header for the sequence

`sed '/^>/ d' Lipotes_vexillifer_noXY_100kb.rounddown.fasta | sed -e '1i\>Baiji' > labeled.fa`
- Prepare mvf file

`python mvftools.py ConvertFasta2MVF --fasta Fasta_files/All_relabeled_combined.fasta --out All_relabeled_combined.mvf`
- Count Dstats patterns

`python /groups/hologenomics/westbury/Software/mvftools/mvftools.py CalcPatternCount --mvf All_relabeled_combined.mvf --out Finless_Harbour_porpoise_Pilotwhale_Orca_Baiji.100kb.input --windowsize 100000 --sample-labels Finless,Harbour_porpoise,Pilotwhale,Orca,Baiji`
- Run Dfoil

`dfoil.py --infile Finless_Harbour_porpoise_Pilotwhale_Orca_Baiji.100kb.input --out Finless_Harbour_porpoise_Pilotwhale_Orca_Baiji.100kb.output`
- Analyse Dfoil for significant patterns of gene flow

`dfoil_analyze.py $file > $file.result`

## Dsuite https://github.com/millanek/Dsuite
- Prepare simulated data for outgroup baiji due to lack of raw reads

`reformat.sh minlength=100000 in=GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna out=GCF_000442215.1_Lipotes_vexillifer_v1_genomic_100kb.fa`

`wgsim -1 150 -2 150 -r 0 -R 0 -X 0 -e 0 -N 100000000 GCF_000442215.1_Lipotes_vexillifer_v1_genomic_100kb.fa reference_1.fastq reference_2.fastq`
- Make vcf

`bcftools mpileup --threads 5 -f GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna -Ou -R Lipotes_100kb.txt -b bamlist.txt | bcftools call -mv -Ov --threads 5 -o All.vcf`
- Calculate the D (ABBA-BABA) and f4-ratio statistics for all possible trios of populations/species

`Dsuite Dtrios -t Tree.txt -o testDtri All_test.vcf Sets.txt`
- Calculate Fbranch

`Dsuite Fbranch Tree.txt testDtri_tree.txt > fbranch.txt`
Plot Fbranch

`dtools.py fbranch.txt Tree.txt`

## Mutation rate https://github.com/ANGSD/angsd

`angsd -doIBS 2 -nthreads 5 -mininddepth 5 -minmapq 30 -minq 30 -uniqueonly 1 -only_proper_pairs 1 -docounts 1 -minInd 9 -makematrix 1 -b bamlist.txt -rf Lipotes_vexillifer_noXY_100kb.txt -ref GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna -out Toothed_whales_IBS` 

## hPSMC https://github.com/jacahill/hPSMC
- Merge 2x fasta consensuses into psmcfa file

`psmcfa_from_2_fastas.py -b10 -m5 species1_autosomes_100kb.fa species2_autosomes_100kb.fa > species1_species2_hPSMC.psmcfa`
- Run psmc

`psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o species1_species2_hPSMC.psmc species1_species2_hPSMC.psmcfa`
- Plot PSMC to get predivergence Ne

`psmc_plot.pl -g 20 -u 7.90E-09 -R -s 10 test species1_species2_hPSMC.psmc`
- run simulations

`hPSMC_quantify_split_time.py -N 20000 -l 500000 -u 5000000 -p 11 -s 11 -o output_prefix`

## Relative dating
`Use scripts found in whaleDatingScripts directory

## Heterozygosity https://github.com/ANGSD/angsd

`angsd -GL 1 -mininddepth 5 -minmapq 30 -minq 30 -uniqueonly 1 -only_proper_pairs 1 -docounts 1 -i bamfile -ref GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna -P 5 -out species -doSaf 1 -anc GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna -rf Lipotes_vexillifer_noXY_100kb.txt -baq 1

`realSFS -nSites 20000000 species.saf.idx -P 10 -tole 1e-8 2> log > species_20Mb.sfs`

## PSMC https://github.com/lh3/psmc

`samtools mpileup -Q 30 -uf reference.fasta bamfile.bam | bcftools call -c - | vcfutils.pl vcf2fq -d 10 | gzip > diploid.fq.gz`

`utils/fq2psmcfa -q20 diploid.fq.gz > diploid.psmcfa`

`utils/splitfa diploid.psmcfa > diploid.split.psmcfa`
- Create PSMC data with whole sequence data

`psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o diploid.psmc diploid.psmcfa`
- Bootstrap replicates (with multiple threads)

`seq 100 | xargs -P 10 -i psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o diploid.round-{}.psmc diploid.split.psmcfa
- Combine the psmc output files

`cat diploid.psmc diploid.round-*.psmc > combined.diploid.psmc`




