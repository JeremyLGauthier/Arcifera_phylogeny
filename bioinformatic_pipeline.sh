#!/bin/bash -       
#title           :bioinformatic_pipeline
#description     :This script contains all command used for the Arcifera phylogeny paper
#author		 : Dr. Jérémy Gauthier
#date            :01.03.2024

#==============================================================================

#RAW READS EVALUATION USING FASTQC
for i in `ls *gz`
	do 
	fastqc "$i"
	done 


#DEMULTIPLEXING
while read a
	do
	name=`echo "$a" | awk '{print $1}'`
	echo "$a" | awk '{print ">"$1"#^"$2}' | tr "#" "\n" > temp_barcodes.fasta
	cutadapt -e 0.17 --max-n 1 --minimum-length 30 -q 10 --no-indels -g file:temp_barcodes.fasta -o demux-"$name"_R1_.fastq.gz -p demux-"$name"_R2_.fastq.gz "$name"_L4_R1_001_*.fastq.gz "$name"_L4_R2_001_*.fastq.gz
	rm temp_barcodes.fasta
	done < list_bar_ind.txt

#CLEANING
for i in `ls demux-*_R1_.fastq.gz`
	do
	cutadapt -u 6 -q 10 -m 30 -a AGATCGGAAGAGC -o clean_"$i" "$i"
	done

for i in `ls demux-*_R2_.fastq.gz`
	do
	cutadapt -u 5 -q 20 -m 30 -a AGATCGGAAGAGC -o clean_"$i" "$i"
	done

#PAIRING
for i in `ls clean_demux-*_R1_.fastq.gz`
		do
		name=`echo $i | sed -e 's/_R1_.fastq.gz//g'`
		gzip -d "$i"
		gzip -d "$name"_R2_.fastq.gz
		fastq_pair "$name"_R1_.fastq "$name"_R2_.fastq
		done

#MAPPING AND CLEANING
for i in `ls *_R1_.fastq.paired.fq`
		do
		sample=`echo $i | sed -e 's/_R1_.fastq.paired.fq//g' -e 's/clean_demux-//g'`
		bwa mem CARAP5_L1.consens_uniq.fa clean_demux-"$sample"_R1_.fastq.paired.fq clean_demux-"$sample"_R2_.fastq.paired.fq > "$sample"_on_ref.sam
		samtools sort "$sample"_on_ref.sam -o temp_1_sorted.bam
		samtools view -bF 4 temp_1_sorted.bam > temp_1_sorted_keep.bam
		source /local/anaconda3/bin/activate /home/jeremy/local/envpicard/
		java -jar -mx128G /home/jeremy/local/envpicard/share/picard-2.20.2-0/picard.jar CreateSequenceDictionary R=CARAP5_L1.consens_uniq.fa O=CARAP5_L1.consens_uniq.dict
		source /local/anaconda3/bin/activate /home/jeremy/local/envpicard/
		java -jar -mx128G /home/jeremy/local/envpicard/share/picard-2.20.2-0/picard.jar AddOrReplaceReadGroups I=temp_1_sorted_keep.bam O=temp_1_sorted_keep_rg.bam ID=["$sample"] RGLB=[id] PL=[pl] PU=[pu] SM=["$sample"]
		samtools index temp_1_sorted_keep_rg.bam
		source /local/anaconda3/bin/activate /home/jeremy/local/envgatk
		GenomeAnalysisTK -T RealignerTargetCreator -I temp_1_sorted_keep_rg.bam -R CARAP5_L1.consens_uniq.fa -o temp.intervals
		java -jar -Xmx128G /home/jeremy/local/envgatk/opt/gatk-3.8/GenomeAnalysisTK.jar -T IndelRealigner -I temp_1_sorted_keep_rg.bam -R CARAP5_L1.consens_uniq.fa -targetIntervals temp.intervals -o "$sample"_sorted_keep_rg_realign.bam
		rm temp*
		rm -r ./results_"$sample"_sorted_keep_pcrdup
		rm CARAP5_L1.consens_uniq.dict
		done
				
#CONSENSUS
for i in `ls *_realign.bam`
	do
	samtools mpileup -A -uf CARAP5_L1.consens.fa "$i" | bcftools call --ploidy 1 -c | bcftools filter -e 'DP<2' | vcfutils.pl vcf2fq > "$i"_consensus_2X.fastq
	done
for i in `ls *_realign.bam`
	do
	samtools mpileup -A -uf CARAP5_L1.consens.fa "$i" | bcftools call --ploidy 1 -c | bcftools filter -e 'DP<3' | vcfutils.pl vcf2fq > "$i"_consensus_3X.fastq
	done
for i in `ls *_realign.bam`
	do
	samtools mpileup -A -uf CARAP5_L1.consens.fa "$i" | bcftools call -c | bcftools filter -e 'DP<6' | vcfutils.pl vcf2fq > "$i"_consensus_6X.fastq
	done

#LOCUS ALIGNMENTS
while read a
	do
	locus=`echo $a | sed -e 's/>//g'` 
	for i in `ls *_2X.fasta`
		do
		sample=`echo $i | sed -e 's/_sorted_keep_rg_realign.bam_consensus_2X.fasta//g'`
		grep -A 1 "$a"$ "$i" | sed 's/>.*/>'"$sample"'/g' >> "$locus"_2X.fasta
		done
	conda activate /data/work/PLATYCARABUS/ANALYSES/SAVE/SANS_PCR_DUP/envmafft
	mafft --auto "$locus"_2X.fasta > "$locus"_2X_align.fasta
	done < temp_list_loci

find . -size 0 -delete

#FILTERS
for i in `ls *_align.fasta`
	do
	name=`echo $i | sed -e 's/.fasta//g'`
	sed -e 's/-/n/g' "$i" > "$i"_2
	seqtk comp "$i"_2 | awk '{print $1"\t"$9/$2}' | awk '{if ($2 < 0.8) print $1}' > temp_list_keep
	fastaselect.pl "$i"_2 temp_list_keep > "$name"_N80.fasta
	done

#CONCAT AND CONVERSION
AMAS.py concat -i *.fasta -f fasta -d dna
sed -e 's/-/N/g' -e 's/?/N/g' concatenated.out > concatenated_2.out
./fasta2relaxedPhylip.pl -f concatenated_2.out -o concatenated_2.phylip

#PHYLOGENY USING IQTREE
#$1 alignment phylip format
#$2 partition file

#PART 1 PARTITION FINDER TO IDENTIFY PARTITIONS
mkdir FOLDER_2
ln -s ../$1 ./FOLDER_2/$1

header="
## ALIGNMENT FILE ##\n
alignment = $1;\n
\n
## BRANCHLENGTHS: linked | unlinked ##\n
branchlengths = linked;\n
\n
## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##\n
models = GTR, GTR+G, GTR+I+G;\n
\n
# MODEL SELECCTION: AIC | AICc | BIC #\n
model_selection = aicc;\n
\n
## DATA BLOCKS: see manual for how to define ##\n
[data_blocks]"

echo -e $header >> ./FOLDER_2/partition_finder.cfg
awk '{print $1";"}' $2 >> ./FOLDER_2/partition_finder.cfg

end="\n
## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans##\n
[schemes]\n
search = rcluster;\n"

echo -e $end >> ./FOLDER_2/partition_finder.cfg

source /local/anaconda3/bin/activate /home/jeremy/local/envforpartitionfiner2
PartitionFinder.py -p 6 -r --min-subset-size=2000 FOLDER_2/

#PART 2 IQ TREE MODEL TO EACH PARTITION
sed -e 's/RaxML/@/g' -e 's/#nexus/@#nexus/g' ./FOLDER_2/analysis/best_scheme.txt | tr '\n' '?' | tr '@' '\n' | grep "#nexus" | tr '?' '\n' | grep -v "charpartition PartitionFinder" > partitionfinder_best_scheme_partitions.nex
conda deactivate

source /local/anaconda3/bin/activate /home/jeremy/local/enviqtree1.6.11/
iqtree -s $1 -m MF -spp partitionfinder_best_scheme_partitions.nex -safe

#PART 3 IQ TREE PHYLOGENY
iqtree -s $1 -spp partitionfinder_best_scheme_partitions.nex.best_scheme.nex -bb 1000 -bnni -alrt 1000 -nt AUTO --runs 100 -allnni -safe -pre concatenated_alignment_output


#ASTRAL PHYLOGENY
for i in `ls *.fasta`
	do
	iqtree -s "$i" -st DNA -m MFP -bb 1000 -bnni -AICc -T AUTO -safe
	done
cat *.treefile > all_treefile
astral-hybrid -i all_treefile -o astral_tree

#SPECIES DELIM
#BPP
bpp --cfile ./platy.bpp.A11_alg0.ctl
          seed =  -1

       seqfile = ./cat_loci_44s.txt
      Imapfile = ./platy.Imap.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

  speciesdelimitation = 1 1 2 0.5  * speciesdelimitation algorithm1 finetune (a m)
          speciestree = 1          * species tree estimation using SPR+SNL

    speciesmodelprior = 1          * 0: uniform labeled histories; 1:uniform rooted trees

  species&tree = 11  HVN HVV HMM CAR CIN PCY PDL PDE PCR PFA PIR
                     3 2 1 2 3 4 3 6 4 7 9 
                    ((HVN,HVV),((HMM,(CAR,CIN)),(PCY,((PDL,PDE),(PCR,(PFA,PIR))))));



       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 221    * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = gamma 2 1000   # Gamma(a,b) for theta
      tauprior = gamma 2 1000   # Gamma(a,b) for root tau

      finetune = 1

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 8000
      sampfreq = 100 
       nsample = 100000

#TR2
for i in `ls *.treefile` ; do pxrr -t "$i" -r -g APCGn1n29,CBX0193,CBX0182,CBX0192,CBX0190 -o "$i"_root.tree ; done
cat *_root.tree > all_treefile_rooted
run_tr2.py -t all_treefile_rooted


#POPULATION GENOMICS PART
#SNP calling
samples=""

for data in `ls *_realign.bam`
	do
	samples=$samples" -I "$data
	done

gatk HaplotypeCaller -R CARAP5_L1.consens.fa -O gatk4.vcf $samples

vcftools --vcf gatk4.vcf --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.4 --recode --out snp_bi_miss04
vcf2structure_gn.sh snp_bi_miss04.recode.vcf

#Loop structure
for i in `seq 1 3`
	do
	for j in `seq 11 15`
		do
		mkdir st"$i"_"$j"_"$1".dir
		nb_sample=`grep -c "." $1 | awk '{print $1/2}'`
		#nb_marker=`awk '{print NF-1}' $1 | sort | uniq -c | awk '{print $2}'`
		nb_marker=`head -n 1 $1 | tr '\t' '\n' | grep "." -c | awk '{print $1-1}'`
		sed -e 's/KN/'"$j"'/g' -e 's/name_output/st'"$i"'_'"$j"'_results/g' -e 's/name_sample/'"$1"'/g' -e 's/nb_sample/'"$nb_sample"'/g' -e 's/nb_marker/'"$nb_marker"'/g' mainparams > ./st"$i"_"$j"_"$1".dir/mainparams
		cp extraparams ./st"$i"_"$j"_"$1".dir/
		cp structure ./st"$i"_"$j"_"$1".dir/
		r1=`shuf -i1-9 -n1`
		r2=`shuf -i1-9 -n1`
		r3=`shuf -i1-9 -n1`
		cp run.sh ./st"$i"_"$j"_"$1".dir/st"$i"_"$j".sh
		sed -e 's/structure/structure\ \-D\ '"$r1"''"$r2"''"$r3"'/g' ./st"$i"_"$j"_"$1".dir/st"$i"_"$j".sh > ./st"$i"_"$j"_"$1".dir/st"$i"_"$j"_2.sh
		cd ./st"$i"_"$j"_"$1".dir/
		chmod 755 ./st"$i"_"$j"_2.sh
		nohup ./st"$i"_"$j"_2.sh &
		cd ..
		done
	done

#DSUITE
Dsuite Dtrios snp_bi_miss04.recode.vcf SETS.txt -t species_tree.nwk
Dsuite Fbranch species_tree.nwk SETS_tree.txt > fbranch.txt
dtools.py fbranch.txt species_tree.nwk

#PHYLOGENY ON SNPs
python vcf2phylip.py -i snp_bi_miss04.recode.vcf
python ascbias.py -p snp_bi_miss04.recode.vcf.phy
raxml-ng --all --msa out.phy --model GTR+ASC_LEWIS --tree pars{10} --bs-trees 100






	
