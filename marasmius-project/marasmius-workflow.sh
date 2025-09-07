#organise data
	mkdir marasmius/ #this is the working directory
	alias marasmius='cd /mnt/sda/proj/Thibault/marasmius/' #shortcut
	mkdir marasmius/data
	cp -rs /mnt/sda/proj/Thibault/raw-data/pr_017/ ./marasmius/data
	#local copy command to transfer to dardel
	scp -r thibault@emb-p-susrv01.emb.su.se:/home/thibault/marasmius/data/pr_017/ ./marasmius/data
	#I rename the strain but keep the original identification numbers of pr_017

#install conda
	https://docs.anaconda.com/free/miniconda/ (download miniconda3, scp on spore)
	bash Miniconda3-latest-Linux-x86_64.sh #installed it in /home/thibault/ but only 18Go, moved it in /mnt/sda/proj/Thibault/miniconda3
	#process to move conda
	conda activate environment #for each environment
	conda env export > environment.yml #export it
	rm -rf ~/miniconda3 #remove it
	bash Miniconda3-latest-Linux-x86_64.sh #reinstall it in the new location
	conda env create -f /home/thibault/environment.yml #install each environment
		conda config --add channels defaults
		conda config --add channels bioconda
		conda config --add channels conda-forge

#assess assemblies in data/analysis
	#busco
		conda create -n busco busco # + activate environement (I won't write it explicitely here)
		cd data/analysis
		for i in $(ls */pr*.fasta); do echo busco -i $i -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 60 --out_path busco_results -q;done >busco.sh
		nohup busco.sh >busco log 2>&1 &
	#quast
		conda create -n quast quast
		mkdir quast_results
		for i in $(ls */pr*.fasta); do quast $i -o quast_results/$i;done > quast.logerr 2>&1 &
#assembly
	mkdir assemblies
	#hifiasm (see other in "not selected")
		conda create -n hifiasm hifiasm
		mkdir assemblies/1_hifiasm
		for i in $(ls ccsreads); do echo hifiasm -o assemblies/1_hifiasm/$i -t 48 ccsreads/$i/$i.ccsreads.fastq.gz \>assemblies/1_hifiasm/$i.log 2\>\&1 ; done | bash &
		cd assemblies/1_hifiasm
		for i in $(ls | grep "p_ctg.gfa"); do awk '/^S/{header=">"$2; for(i=4; i<=NF; i++) {header=header" "$i}; print header; printf "%s", $3 | "fold -w 80"; close("fold -w 80"); print ""}' $i > $i\.fa;done
		#exact same results as in analysis, I take those as they are the best I got
		#busco
			mkdir busco_results
			for i in $(ls | grep gfa.fa | grep -v ".fai"); do echo busco -i $i -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 60 --out_path busco_results \>busco_results/$i.log 2\>\&1; done > busco.sh
			nohup bash busco.sh 2>busco.log 2>&1 & [2] 2056537
			#for i in $(ls | grep BUSCO); do sed -n 9p $i/*.txt; done > busco_results.txt #formatting in calc
		#quast
			mkdir quast_results
			for i in $(ls *.fa); do quast $i -o quast_results/$i;done > quast.logerr 2>&1 &
		#I take the best haplotype or whole primary assembly for each strain, based on the busco results
		cp baemy1.bp.hap1.p_ctg.gfa.fa baemy1.hifiasm.fa
		cp colpe1.bp.p_ctg.gfa.fa colpe1.hifiasm.fa
		cp gymaq1.bp.hap1.p_ctg.gfa.fa gymaq1.hifiasm.fa
		cp marsi1.bp.p_ctg.gfa.fa marsi1.hifiasm.fa
		cp mycal1.bp.hap1.p_ctg.gfa.fa mycal1.hifiasm.fa
		cp mycal2.bp.hap1.p_ctg.gfa.fa mycal2.hifiasm.fa
		cp mycsc1.bp.hap1.p_ctg.gfa.fa mycsc1.hifiasm.fa
		cp rhodo1.bp.hap1.p_ctg.gfa.fa rhodo1.hifiasm.fa
	#illumina assemblies
		#pooling runs
			for i in $(ls); do gunzip $i/*;done &
			cat Sample_WE-3693-M1-Gymer/WE-3693-M1-Gymer_S112_L003_R1_001.fastq Sample_WE-3693-M1-Gymer_1/WE-3693-M1-Gymer_S113_L003_R1_001.fastq > WE-3693-M1-Gymer_S112-S113_L003_R1_001.fastq
			cat Sample_WE-3693-M1-Gymer/WE-3693-M1-Gymer_S112_L003_R2_001.fastq Sample_WE-3693-M1-Gymer_1/WE-3693-M1-Gymer_S113_L003_R2_001.fastq > WE-3693-M1-Gymer_S112-S113_L003_R2_001.fastq
			cat Sample_WE-3693-M21-Marro/WE-3693-M21-Marro_S114_L003_R1_001.fastq Sample_WE-3693-M21-Marro_1/WE-3693-M21-Marro_S115_L003_R1_001.fastq > WE-3693-M21-Marro_S114-S115_L003_R1_001.fastq
			cat Sample_WE-3693-M21-Marro/WE-3693-M21-Marro_S114_L003_R2_001.fastq Sample_WE-3693-M21-Marro_1/WE-3693-M21-Marro_S115_L003_R2_001.fastq > WE-3693-M21-Marro_S114-S115_L003_R2_001.fastq
			gzip -9 WE-3693-M1-Gymer_S112-S113_L003_R1_001.fastq &
			gzip -9 WE-3693-M1-Gymer_S112-S113_L003_R2_001.fastq &
			gzip -9 WE-3693-M21-Marro_S114-S115_L003_R1_001.fastq &
			gzip -9 WE-3693-M21-Marro_S114-S115_L003_R2_001.fastq &
		#fastqc
		conda create -n fastqc fastqc
		cd data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7
		nohup fastqc -f fastq -noextract -t 48 WE-3693-M1-Gymer_S112-S113_L003_R1_001.fastq.gz WE-3693-M1-Gymer_S112-S113_L003_R2_001.fastq.gz > fastqc_m1.logerr 2>&1 &
		nohup fastqc -f fastq -noextract -t 48 WE-3693-M21-Marro_S114-S115_L003_R1_001.fastq.gz WE-3693-M21-Marro_S114-S115_L003_R2_001.fastq.gz > fastqc_M21.logerr 2>&1 &
		#locally
		scp thibault@emb-p-susrv01.emb.su.se:/mnt/sda/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M1-Gymer_S112-S113_L003_R1_001_fastqc.html .
		scp thibault@emb-p-susrv01.emb.su.se:/mnt/sda/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M1-Gymer_S112-S113_L003_R2_001_fastqc.html .
		scp thibault@emb-p-susrv01.emb.su.se:/mnt/sda/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M21-Marro_S114-S115_L003_R1_001_fastqc.html .
		scp thibault@emb-p-susrv01.emb.su.se:/mnt/sda/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M21-Marro_S114-S115_L003_R2_001_fastqc.html .
		#MaSuRCA
			conda create -n masurca masurca
			mkdir assemblies/masurca
			mkdir assemblies/masurca/gymer
			cd assemblies/masurca/gymer
			nohup masurca -i /mnt/sda/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M1-Gymer_S112-S113_L003_R1_001.fastq.gz, /mnt/sda/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M1-Gymer_S112-S113_L003_R2_001.fastq.gz -t 16 >/mnt/sda/proj/Thibault/marasmius/assemblies/masurca/gymer/masurca.logerr 2>&1 &
			mkdir assemblies/masurca/marro
			nohup masurca -i /mnt/sda/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M21-Marro_S114-S115_L003_R1_001.fastq.gz, /mnt/sda/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M21-Marro_S114-S115_L003_R2_001.fastq.gz -t 16 >/mnt/sda/proj/Thibault/marasmius/assemblies/masurca/marro/masurca.logerr 2>&1 &
			#busco
				cd masurca/gymer/ #same for marro
				mkdir busco_results/
				busco -i CA/primary.genome.scf.fasta -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 60 --out_path busco_results -q > busco_results/busco.log 2>&1 &
			#quast
				cd masuca/gymer
				mkdir quast_results/
				quast CA/primary.genome.scf.fasta -o quast_results/ >quast_results/quast.out & #same for marro
			#I keep those as they are the best I got
	#quast
		conda create -n quast quast
		cd marasmius/data/analysis
		mkdir quast_results
		for i in $(ls); do echo quast $i/* -o quast_results/$i; done | bash
		cd quast_results/
		for i in $(ls); do tail -n 1 $i/transposed_report.tsv;done > transposed_report.tsv #paste header and transpose via calc
	#assembly decontamination blobtools
conda create -n blobtools blobtools
	#mapping
		conda create -n samtools samtools
		conda create -n minimap2 minimap2
		mkdir 2.2_blob-hap1/
		mkdir 2.2_blob-hap1/mapping
		nohup minimap2 -a -t 60 -x map-hifi ../1_hifiasm/baemy1.hifiasm.fa ../../data/ccsreads/baemy1_004/pr_017_004.ccsreads.fastq.gz -o mapping/baemy1.sam > mapping/baemy1.stderr 2>&1 & #same for others
		#for illumina (same for marro)
			cd assemblies/masurca/gymer/CA
			conda create -n bowtie2 bowtie2
			bowtie2-build primary.genome.scf.fasta gymer >bowtie2-build.log 2>&1
			cd marasmius
			nohup bowtie2-align-s -x assemblies/masurca/gymer/CA/gymer -1 data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M1-Gymer_S112-S113_L003_R1_001.fastq.gz -2 data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M1-Gymer_S112-S113_L003_R2_001.fastq.gz -S mapping/gymer.sam -p 48 > mapping/gymer-bowtie2.log 2>&1 &
		cd mapping
		for i in $(ls *.sam); do out=$(echo $i | sed "s/sam/_sorted.bam/"); samtools sort -@ 12 -o $out $i;done >sorttobam.log 2>&1
		rm *.sam
		for i in $(ls *_sorted.bam); do samtools index $i;done
	#hit
		#split the genomes
			cd assemblies/1_hifiasm
			conda create -n seqkit seqkit
			seqkit split -s 1 baemy1.hifiasm.fa -O ../2.2_blob-hap1/hit/baemy1/ #same for others
			seqkit split -p 96 /assemblies/masurca/gymer/CA/primary.genome.scf.fasta # for illumina
		#diamond blastx
			cd 2_blob/
			mkdir hit/
			for i in $(ls *.fa); do echo diamond blastx -p 18 -d /mnt/sdb/blast_databases/nr.dmnd -q $i -f 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore -o $i\_blastx.tsv \>  $i\_blastx.log 2\>\&1;done > baemy1.sh
			nohup bash baemy1.sh > baemy1.log 2>&1 & #same for others, I ran in ultra-sensitive mode the ones that show no results
		#merge and check results
			cat *tsv> baemy1.hifiasm_blastx.out #same for others
			for i in $(ls *tsv); do echo $i;cut -f2 $i;done #check that the taxid is here
			gawk -i inplace -v INPLACE_SUFFIX=.bak -v FS=$'\t' -v OFS=$'\t' '{if ($2) print $0;}' hit/colpe1/colpe1.bp.hap1.p_ctg.gfa_blastx.tsv #remove lines with empty taxid
	#blob
		wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
		gunzip taxdump.tar.gz
		tar -xf taxdump.tar
		conda create -n blobtools blobtools
		blobtools create -i ../1_hifiasm/baemy1.hifiasm.fa -b mapping/baemy1_sorted.bam -t hit/baemy1/baemy1.hifiasm_blastx.tsv --names names.dmp -o baemy1 #nothing to remove
		blobtools plot -i baemy1.blobDB.json #same for others
		blobtools view -i baemy.blobDB.json #same for others
		blobtools create -i ../1_hifiasm/gymaq1.bp.hap1.p_ctg.gfa.fa -b mapping/gymaq1_sorted.bam -t hit/gymaq1/gymaq1.hifiasm_blastx.tsv --names names.dmp -o gymaq1 #two to remove
		blobtools create -i ../1_hifiasm/mycal1.bp.hap1.p_ctg.gfa.fa -b mapping/mycal1_sorted.bam -t hit/mycal1/mycal1.hifiasm_blastx.tsv --names names.dmp -o mycal1 #many to remove
		blobtools create -i ../1_hifiasm/mycal2.bp.hap1.p_ctg.gfa.fa -b mapping/mycal2_sorted.bam -t hit/mycal2/mycal2.hifiasm_blastx.tsv --names names.dmp -o mycal2 #one to remove
		blobtools create -i ../1_hifiasm/mycsc1.bp.hap1.p_ctg.gfa.fa -b mapping/mycsc1_sorted.bam -t hit/mycsc1/mycsc1.hifiasm_blastx.tsv --names names.dmp -o mycsc1 #one to remove
		blobtools create -i ../1_hifiasm/rhodo1.bp.hap1.p_ctg.gfa.fa -b mapping/rhodo1_sorted.bam -t hit/rhodo1/rhodo1.hifiasm_blastx.tsv --names names.dmp -o rhodo1 #two to remove
		#download plots (local, same for others)
	#cleaning
		#baemy1: nothing to clean
		awk -vFS=$'\t' -v OFS=$'\t' 'NR>11{if($6 == "Basidiomycota" || $6 == "Ascomycota"){print $1}}' gymaq1.blobDB.table.txt > gymaq1.clean
		awk -vFS=$'\t' -v OFS=$'\t' 'NR>11{if($6 == "Basidiomycota" || $6 == "Ascomycota"){print $1}}' marsi1.blobDB.table.txt > marsi1.clean
		awk -vFS=$'\t' -v OFS=$'\t' 'NR>11{if($6 == "Basidiomycota" || $6 == "Ascomycota"){print $1}}' mycal1.blobDB.table.txt > mycal1.clean #mycal1
		awk -vFS=$'\t' -v OFS=$'\t' 'NR>11{if($6 == "Basidiomycota" || $6 == "Ascomycota" || $6 == "no-hit"){print $1}}' mycal2.blobDB.table.txt > mycal2.clean
		awk -vFS=$'\t' -v OFS=$'\t' 'NR>11{if($6 == "Basidiomycota" || $6 == "Ascomycota"){print $1}}' mycsc1.blobDB.table.txt > mycsc1.clean
		awk -vFS=$'\t' -v OFS=$'\t' 'NR>11{if($6 == "Basidiomycota" || $6 == "no-hit"){print $1}}' rhodo1.blobDB.table.txt > rhodo1.clean
		awk 'NR>11{if($6 =="Basidiomycota" || $6 =="Ascomycota" || ($6 =="no-hit" && $3 >= 0.28))print $1}' gymer.blob4.blobDB.table.txt > gymer.clean #in 2_blob/
		awk 'NR>11{if($6 =="Basidiomycota" || $6 =="Ascomycota" || $6 =="no-hit")print $1}' marro.blob4.blobDB.table.txt > marro.clean #in 2_blob/
		cd assemblies/1_hifiasm
		blobtools seqfilter -i gymaq1.hifiasm.fa -l ../2.2_blob-hap1/gymaq1.clean
		blobtools seqfilter -i marsi1.hifiasm.fa -l ../2.2_blob-hap1/marsi1.clean
		blobtools seqfilter -i mycal1.hifiasm.fa -l ../2.2_blob-hap1/mycal1.clean
		blobtools seqfilter -i mycal2.hifiasm.fa -l ../2.2_blob-hap1/mycal2.clean
		blobtools seqfilter -i mycsc1.hifiasm.fa -l ../2.2_blob-hap1/mycsc1.clean
		blobtools seqfilter -i rhodo1.hifiasm.fa -l ../2.2_blob-hap1/rhodo1.clean
		cd assemblies/1_masurca/gymer/CA/ #same for marro
		blobtools seqfilter -i primary.genome.scf.fasta -l ../../../../blob/gymer.clean
	#busco
		cd 1_masurca/gymer/
		busco -i CA/primary.genome.scf.filtered.fna -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 60 --out_path busco_results -q > busco_results/busco.filtered.log 2>&1 &
		cd 1_masurca/marro/
		busco -i CA/primary.genome.scf.filtered.fna -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 60 --out_path busco_results -q > busco_results/busco.filtered.log 2>&1 & 
		cd 1_hifiasm
		for i in $(ls *filtered.fna); do echo busco -i $i -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 60 --out_path busco_results -q; done | bash > busco_results/busco.filtered.log 2>&1 & 
	for i in $(ls | grep filtered); do mv $i $(echo $i | sed 's/filtered.fna/blob.fa/');done #rename the output in blob.fa
	cp baemy1.hifiasm.fa baemy1.hifiasm.blob.fa #just for consistency
	#purge_dups
		cd assemblies/
		mkdir 3_purge_dups
		conda create -n purge_dups
		for i in $(ls /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/*blob.fa);do strain=$(echo $i | sed 's#.*/\([^/.]*\)\..*#\1#'); echo pd_config.py -n $strain\.config.json -l /mnt/sdb/proj/Thibault/marasmius/assemblies/3_purge_dups/$strain $i /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/$strain*/*fastq.gz;done>pd_config.sh
		sed -i 's/"lineage": *"[^"]*"/"lineage": "agaricales"/' *.json #change whatever lineage to agaricales
		run_purge_dups.py /mnt/sdb/proj/Thibault/marasmius/assemblies/purge_dups/baemy1.config.json /mnt/sdb/proj/Thibault/miniconda3/envs/purge_dups/bin baemy1 > baemy1.purge_dups.log 2>&1 & #brings an error with plit_fa, the specific command works alone though, I switch to the full pipeline guide
		#script for baemy (same for others)
			eval "$(conda shell.bash hook)"
			conda activate purge_dups
			minimap2 -xasm20 /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/baemy1_004/pr_017_004.ccsreads.fastq.gz baemy1.hifiasm.blob.fa -t 40 | gzip -c - > baemy1.paf.gz
			pbcstat baemy1.paf.gz >pbcstat.log 2>&1
			calcuts -d 1 PB.stat > cutoffs 2> calcuts.log #Should be considered as haploid
			split_fa baemy1.hifiasm.blob.fa > baemy1.hifiasm.blob.fa.split
			minimap2 -xasm5 -DP baemy1.hifiasm.blob.fa.split baemy1.hifiasm.blob.fa.split -t 40 | gzip -c - > baemy1.hifiasm.blob.fa.split.self.paf.gz
			purge_dups -2 -T cutoffs -c PB.base.cov baemy1.paf.gz > dups.bed 2>purge_dups.log #should be on the split.self.paf but only work like this?
			get_seqs -e dups.bed baemy1.hifiasm.blob.fa #does not make any difference without the "e" option
			mkdir busco_results
			busco -i purged.fa -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 20 --out_path busco_results -q > busco_results/purged.log 2>&1 
			conda deactivate
			conda deactivate
		nohup bash baemy1_purge-dups.sh 2>&1 & #run analog script for others
		#purge_dups decreases a bit the busco results, but mtLink overcompensates it
	#scaffolding with ntLink
		cd assemblies
		mkdir ntlink
		mkdir ntlink/baemy1 #same for others
		cp ../3_purge_dups/baemy1/purged.fa ntlink/baemy1/baemy1.hifiasm.blob.purged.fa #same for others
		conda create -n ntlink ntlink
		for i in $(ls */ -d); do strain=$(echo $i|sed 's#/##');echo ntLink  scaffold gap_fill target=$strain/$strain.hifiasm.blob.purged.fa reads= ../../data/ccsreads/$strain*/*.fastq.gz t=8 \>$strain/$strain.ntlink.log 2\>\&1 | sed 's#reads\= #reads\=#';done >ntlink.sh
		for i in $(ls */*fill.fa); do echo busco -i $i -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 60 --out_path busco_results -q;done| bash > busco.log 2>&1 &
#available genomes                                      
	#retrieve (same for others)
	mkdir genomes/ #+ make a subdirectory for each of them
	wget https://github.com/glarue/jgi-query/blob/main/jgi-query.py                    
	awk '{if(!/>/){print toupper($0)}else{print $0}}' #remove soft-masking when necessary
	cp GCA_000143185.2_Schco3_genomic.fna schco1.raw.fa #rename them
#sorting and masking (same for others, in "genomes" and in "assemblies/funannotate")
	#install funannotate (see in "annotation")
	cd assemblies
	mkdir funannotate
	cp 4_ntlink/baemy1/baemy1.hifiasm.blob.purged.fa.k32.w100.z1000.ntLink.scaffolds.gap_fill.fa 5_funannotate/baemy1.hifiasm.blob.purged.ntlink.fa #same for others
	cp 1_masurca/gymer/CA/primary.genome.scf.filtered.fna 5_funannotate/gymer1.masurca.blob.fa #same for marro
	cd funannotate/
	funannotate sort -i baemy1.hifiasm.blob.purged.ntlink.fa -o baemy1.hifiasm.blob.purged.ntlink.sort.fa -b baemy1_s --minlen 0 > baemy1.sort.log 2>&1 & #sort
	nohup funannotate mask -i baemy1.hifiasm.blob.purged.ntlink.sort.fa -m repeatmasker -m repeatmodeler -s agaricales -l /mnt/sdb/proj/Thibault/marasmius/db/repbase_all_2022-11.fasta --cpus 10 -o baemy1.hifiasm.blob.purged.ntlink.sort.masked.fa > baemy1.mask.log 2>&1 & #mask
	mkdir busco_results/
	for i in $(ls *masked.fa); do echo busco -i $i -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 60 --out_path busco_results \>busco_results/$i.busco.log;done> busco.sh 
#cleaned
	cd assemblies/
	mkdir cleaned/
	cp 5_funannotate/baemy1.hifiasm.blob.purged.ntlink.sort.masked.fa 6_cleaned/baemy1.fa #same for others
	cp 5_funannotate/gymer1.masurca.blob.sort.masked.fa 6_cleaned/gymer1.fa #same for marro
	cd genomes/
	mkdir cleaned/
	cp collybiopsis-confluens/colco1.raw.sort.masked.fa cleaned/colco1.fa #same for others
	#quast
	mkdir quast_results
	for i in $(ls *fa); do quast $i -o quast_results/$i;done > quast.logerr 2>&1 &
	cd quast_results/
	for i in $(ls); do tail -n 1 $i/transposed_report.tsv;done > transposed_report.tsv
#iqtree
	#prep input
		cd marasmius
		mkdir iqtree
		mkdir iqtree/
		cd iqtree/
		mkdir busco/
		for i in $(ls /mnt/sdb/proj/Thibault/marasmius/genomes/cleaned/*fa); do echo busco -i $i -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 --augustus --augustus_species coprinus_cinereus -c 10 --out_path ./ -q \>$(echo $i | sed 's#.*/\([^/.]*\).*#\1\.log#') 2\>\&1;done>busco.genomes.sh
		for i in $(ls /mnt/sdb/proj/Thibault/marasmius/assemblies/6_cleaned/*fa); do echo busco -i $i -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 --augustus --augustus_species coprinus_cinereus -c 10 --out_path ./ -q \>$(echo $i | sed 's#.*/\([^/.]*\).*#\1\.log#') 2\>\&1;done>busco.assemblies.sh
		cd ../
		conda create -n buscophylo -c bioconda busco_phylogenomics
		nohup busco_phylogenomics.py -i busco/ -o input_nt_90 --nt -t 40 -psc 90 > buscophylo_nt_90.log 2>&1 & 
	#iqtree
		nohup iqtree -s nt_90/supermatrix/SUPERMATRIX.phylip --alrt 1000 -bb 1000 -nt 40 >iqtree_nt90.log 2>iqtree_nt90.err & #alrt stands for approximate likelihood ratio test, bb for ultra-fastbootstrap
		cd input_nt_90/
		iqtree -t supermatrix/SUPERMATRIX.phylip.treefile --gcf gene_trees_single_copy/ALL.tree -s supermatrix/SUPERMATRIX.phylip -p supermatrix/trimmed_alignments/ -T 20 >gcf.log 2>&1 & #add the gene concordance factor
		nohup iqtree -t supermatrix/SUPERMATRIX.phylip.treefile.cf.tree -s supermatrix/SUPERMATRIX.phylip --scfl 100 -T 20 > scf.log 2>&1 & #add the site-concordance factor
		#no alrt, boostrap, gCF nor sCF for marro and sister, checking why
		iqtree -t supermatrix/SUPERMATRIX.phylip.treefile -s supermatrix/SUPERMATRIX.phylip --gcf gene_trees_single_copy/ALL.tree --scf 100 --prefix check_marro -T 50 >check_marro.log 2>&1 &
#astral
	mkdir astral
	conda create -n astral -c bioconda astral-tree 
	nohup astral -i ../iqtree/nt_90/gene_trees_single_copy/ALL.tree -o nt90-astral.tree 2>astral.log &
#annotation
	#install funnanotate
		https://funannotate.readthedocs.io/en/latest/install.html
		#conda create -n funannotate "python>=3.6,<3.9" funannotate
		#pip3 install distro
		#conda activate funannotate
		mamba create -n funannotate2 funannotate #using conda causes a lot of trouble, python<3.9 is deprecated, missing distro, repeatmasker, repeatmodeler, iprscan
	echo $fasta
		conda install repeatmasker repeatmodeler
	#database
		funannotate setup -i all -b basidiomycota -d ./db
		echo export 'funannotate_db=/mnt/sdb/proj/thibault/marasmius/annotation/funannotate/db' >>~/.bash_profile
		source ~/.bash_profile
		conda activate funannotate   
	#check
		funannotate check --show-versions
		#error: emapper.py not installed
		#error: gmes_petap.pl not installed
		#error: signalp not installed
	#install genemark
		http://topaz.gatech.edu/genemark/license_download.cgi (download linux 64 kernel 3.10 - 5 plus keys)
		scp downloads/gm* thibault@emb-p-susrv01.emb.su.se:/mnt/sdb/proj/thibault/ #locally
		cd /mnt/sdb/proj/thibault/
		tar -xzf gmes_linux_64_4.tar.gz
		gunzip gm_key_64.gz
		cp gm_key_64 ~/.gm_key
		export perl5lib=$conda_prefix/bin/perl #not permanent! (see below in the test)
		cd gmes_linux_64_4/
		bash check_install.bash #caution! it printed completed but funannotate check still showed it missing (see below in the test)
		gmes_petap.pl --es --fungus --cpus 24 --sequence test.softmasked.fa
		#-bash: /mnt/sdb/proj/thibault/gmes_linux_64_4/gmes_petap.pl: /mnt/sdb/proj/thibault/miniconda3/envs/funannotate2/bin/perl: bad interpreter: No such file or directory
		#perl path pointed to an old conda environment
		sed -i 's/funannotate2/funannotate/g' ./*.pl
		echo 'export genemark_path=/mnt/sdb/proj/thibault/gmes_linux_64_4' >>~/.bash_profile #if not found in path
		echo 'export path=$path:/mnt/sdb/proj/thibault/gmes_linux_64_4' >>~/.bash_profile #I put it in the path to be sure
	#signalp
		#https://services.healthtech.dtu.dk/services/signalp-4.1/ #linux version
		scp downloads/signalp-5.0b.linux.tar.gz thibault@emb-p-susrv01.emb.su.se:/mnt/sdb/proj/Thibault/ #locally
		cd /mnt/sdb/proj/thibault/
		bin/signalp -fasta test/euk10.fsa -org euk -format short -prefix euk_10_short #sort test
		#2025/03/21 13:52:53 asset: open /mnt/sdb/proj/thibault/signalp-5.0b/bin/bin/signalp: no such file or directory  #seriously? copying the link
		mkdir bin/bin
		cp bin/signalp bin/bin/
		#test completed
		echo 'export path=$path:/mnt/sdb/proj/thibault/signalp-5.0b/bin' >> ~/.bash_profile
		conda deactivate
		source ~/.bash_profile
	#eggnogg
		conda install -c bioconda eggnog-mapper
		cd marasmius/annotation/funannotate
		download_eggnog_data.py --data_dir ./db >eggnogg.db.log 2>&1 &
		emapper.py #diamond database /mnt/sdb/proj/thibault/miniconda3/envs/funannotate2/lib/python3.9/site-packages/data/eggnog_proteins.dmnd not present. Use download_eggnog_database.py to fetch it
		#the db is not in the expected folder, creating soft links
		mkdir /mnt/sdb/proj/thibault/miniconda3/envs/funannotate/lib/python3.9/site-packages/data
		ln -s /mnt/sdb/proj/thibault/marasmius/annotation/funannotate/db/eggnog.db /mnt/sdb/proj/Thibault/miniconda3/envs/funannotate/lib/python3.9/site-packages/data/
		ln -s /mnt/sdb/proj/thibault/marasmius/annotation/funannotate/db/eggnog_proteins.dmnd /mnt/sdb/proj/Thibault/miniconda3/envs/funannotate/lib/python3.9/site-packages/data/
		ln -s /mnt/sdb/proj/thibault/marasmius/annotation/funannotate/db/eggnog.taxa.db /mnt/sdb/proj/Thibault/miniconda3/envs/funannotate/lib/python3.9/site-packages/data/
		ln -s /mnt/sdb/proj/thibault/marasmius/annotation/funannotate/db/eggnog.taxa.db.traverse.pkl /mnt/sdb/proj/Thibault/miniconda3/envs/funannotate/lib/python3.9/site-packages/data/
		cd test/
		emapper.py -m diamond -i test-rna_seq_fead2739-ec13-4cd5-9c42-05516ab11a4c/rna-seq/predict_misc/proteins.combined.fa --data_dir ../db/ --output_dir test-eggnog -o test --override --cpu 40 >test-eggnog.log 2>&1 & #ok
	#test
		mkdir test
		funannotate test -t all
		#running `funannotate clean` unit testing: minimap2 mediated assembly duplications
		#success: `funannotate clean` test complete.
		#running `funannotate mask` unit testing: repeatmodeler --> repeatmasker
		#success: `funannotate mask` test complete.
		#running `funannotate predict` unit testing
		#cmd error: /mnt/sdb/proj/thibault/gmes_linux_64_4/gmes_petap.pl --es --max_intron 3000 --soft_mask 2000 --cores 2 --sequence genome.query.fasta --fungus
		#[mar 21 01:19 pm]: can't locate findbin.pm in @inc (you may need to install the findbin module) (@INC contains: /usr/local/lib64/perl5/5.32 /usr/local/share/perl5/5.32 /usr/lib64/perl5/vendor_perl /usr/share/perl5/vendor_perl /usr/lib64/perl5 /usr/share/perl5) at /mnt/sdb/proj/Thibault/gmes_linux_64_4/gmes_petap.pl line 79.
		#begin failed--compilation aborted at /mnt/sdb/proj/thibault/gmes_linux_64_4/gmes_petap.pl line 79.
		./change_path_in_perl_scripts.pl /mnt/sdb/proj/thibault/miniconda3/envs/funannotate/bin/perl
		funannotate test -t all -cpu 80 >funannotate.test.log2 #everything works except 'compare' due to iqtree missing dependency, no need fot it though
	#gene prediction
		mkdir baemy1/ #same for others
		for i in $(ls ../../assemblies/6_cleaned/*.fa); do echo funannotate predict -i $i -o $(echo $i | sed -E 's#.*/([^/]+)\.fa#\1/#') -s $(echo $i | sed -E 's#.*/([^/]+)\.fa#"\1"#') --busco_db basidiomycota --ploidy 2 --cpus 80 \>$(echo $i | sed -E 's#.*/([^/]+)\.fa#\1.predict.log#') 2\>\&1;done>predict.assemblies.sh
		for i in $(ls ../../genomes/cleaned/*.fa); do echo funannotate predict -i $i -o $(echo $i | sed -E 's#.*/([^/]+)\.fa#\1/#') -s $(echo $i | sed -E 's#.*/([^/]+)\.fa#"\1"#') --busco_db basidiomycota --ploidy 2 --cpus 80 \>$(echo $i | sed -E 's#.*/([^/]+)\.fa#\1.predict.log#') 2\>\&1;done>predict.genomes.sh
	#fonctionnal annotation
		#interproscan (+eggnog and signalp)
			cd funannotate/db
			wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.59-91.0/interproscan-5.59-91.0-64-bit.tar.gz.md5
			wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.59-91.0/interproscan-5.59-91.0-64-bit.tar.gz 
			tar xvzf interproscan-5.59-91.0-64-bit.tar.gz
			cd interproscan-5.73-104.0/
			python3 setup.py -f interproscan.properties
			for i in $(ls -d */ | grep -Ev '^(test|db)/?$'); do echo funannotate iprscan -i $i -m local -c 20 --iprscan_path /mnt/sdb/proj/Thibault/interproscan-5.73-104.0/ \>$(echo $i | sed 's#/#.iprscan.log#') 2\>\&1;done>iprscan.sh
		#eggnog
			for i in $(ls -d */ | grep -Ev '^(test|db)/?$'); do echo funannotate annotate -i $i -o $i --busco_db basidiomycota --cpus 40 \>$(echo $i | sed 's#/#.annotate.log#') 2\>\&1;done>annotate.sh
		#antismash
			conda create -n antismash antismash
			cd annotation/funannotate
			download-antismash-databases
			pip install bcbio-gff==0.6.9
			pip install bcbio-gff==0.7.0 #instead of 0.7.1
			pip install moods-python==1.9.4.1 #instead of 1.9.4.2
			pip install biopython==1.78 #instead of 1.79
			pip install brawn==1.0.1 #instead of 1.0.2
			pip install jinja2==3.1.2 #instead of 3.1.6
			pip install joblib==1.3.2 #instead of 1.4.2
			pip install jsonschema==4.11.0 #instead of 4.23.0
			pip install markupsafe==2.1.3 #instead of 3.0.2
			pip install matplotlib==3.8.1 #instead of 3.10.1
			pip install numpy==1.26.2 #instead of 2.2.4
			pip install scikit-learn==1.3.2 #instead of 1.6.1
			pip install requires scipy==1.11.3 #instead of 1.15.2
			for i in $(ls -d */ | grep -Ev '^(test|db)/?$'); do echo antismash -i "$i"predict_results/$(echo $i | sed 's#/#.gbk#') -t fungi -c 5 --genefinding-tool glimmerhmm --genefinding-gff3 "$i"predict_results/$(echo $i | sed 's#/#.gff3#') --output-dir $i/antismash --output-basename $(echo $i | sed 's#/##') \>$(echo $i | sed 's#/#.antismash.log#') 2\>\&1;done>antismash.sh
		#phobius
			#downloaded phobius + the sets at https://phobius.sbc.su.se/data.html
			cd /mnt/sdb/proj/Thibault/
			tar -xvzf phobius101_linux.tgz #tar -xvzf set.tar.gz in phobius directory
			cd /mnt/sdb/proj/Thibault/marasmius/annotation/funannotate/baemy1
			conda create -n perl-findbin bioconda::perl-findbin #adding missing dependancy
			mkdir phobius
			seqkit split --by-id --by-id-prefix '' annotate_misc/genome.proteins.fasta -O phobius/ #split the proteins
			perl phobius.pl ../marasmius/annotation/funannotate/baemy1/phobius/FUN_000001-T1.fa #test
			#test can't find the fasta, took some time to figure out, a forum suggested to change decodeanhmm to 64.bit
			vim phobius.pl #line 25
			#test throws "Use of uninitialized value $predstr in concatenation (.) or string at phobius.pl line 216, <SOUT> line 15."
			vim phobius.pl #setting $predstr=""; line 181
			#test with -short option throws the same error, setting setting $predstr=""; line 243
			for i in $(ls *fa); do perl /mnt/sdb/proj/Thibault/phobius/phobius.pl -short $i > $(echo $i | sed 's/fasta/phobius/'); done > phobius.log 2>phobius.err & #throws citation to .err
			#paste everything into a script phobius.sh for the other samples:
				eval "$(conda shell.bash hook)"
				mkdir phobius
				conda activate seqkit
				seqkit split --by-id --by-id-prefix '' annotate_misc/genome.proteins.fasta -O phobius/
				conda deactivate
				#splitMfasta.pl --minsize=0 --outputpath=./phobius/ annotate_misc/genome.proteins.fasta
				#for i in $(ls phobius/); do mv phobius/$i phobius/$(head -n1 $i | cut -f1 -d' ' | sed 's/>//' | sed 's/$/.fa/'); done
				conda activate perl-findbin
				for i in $(ls phobius/*fasta); do perl /mnt/sdb/proj/Thibault/phobius/phobius.pl -short $i > $(echo $i | sed 's/fasta/phobius/'); done
				conda deactivate
			cp phobius.sh
			bash phobius.sh > phobius.log 2>phobius.err & #same for others
			tail -q -n1 *phobius > baemy1.phobius #concatenate results, paste header manually (should have put it in the previous script)
		#combine
			for i in $(ls -d */ | grep -v test |grep -v db |grep -v old|sed 's#/##');do echo funannotate annotate -i $i/ -o $i/ --antismash $i/antismash/$i.gbk --phobius $i/phobius/$i.phobius --cpus 20 \>$i/$i.combine.log 2\>\&1;done >combine.sh
			#warning: the automatic logs (funannotate-annotate.xxx + genome.fixedproducts) are printed in the current working directory, because of the loop i.e. in the parent directory funannotate/, 
		#everything has been put together in a one scrit "funannotate.all.assemblies.sh" containing all the steps
		#curating
		head */annotate_results/Gene2Products.need-curating.txt > curating.txt #don't know what to do with this, I stop
	#earlgrey
		conda create -n earlgrey earlgrey
		conda activate earlgrey
		conda install bioconda::perl-findbin #needs perl-findbin
		wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.16.h5.gz
		wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.16.h5.gz.md5
		mv dfam39_full.16.h5 miniconda3/envs/earlgrey/share/RepeatMasker/Libraries/famdb/
		cd miniconda3/envs/earlgrey/share/RepeatMasker
		perl ./configure #ok for trf path, choose rmblast (hmmer + dfam does not work)
		sed -i '1s|^#!.*perl.*|#!/usr/bin/env perl|' /mnt/sdb/proj/Thibault/miniconda3/envs/earlgrey/share/earlgrey-6.0.0-0/scripts/LTR_FINDER_parallel #after a test LTR_FINDER takes perl fro m the system, change the shebang to use the conda version with perl-findbin
		mkdir earlgrey/
		for i in $(ls ../../assemblies/6_cleaned/*fa); do echo earlGrey -g $i -s $(echo $i | sed -E 's#.*/(.*)\.fa#\1#') -o $(echo $i | sed -E 's#.*/(.*)\.fa#./\1/#') -t 50 -r agaricales -l /mnt/sdb/proj/Thibault/marasmius/db/repbase_all_2022-11.fasta -c yes -e yes \>$(echo $i | sed -E 's#.*/(.*)\.fa#\1.earlgrey.log#') 2\>\&1;done >earlgrey.assemblies.sh
		for i in $(ls ../../genomes/cleaned/*fa); do echo earlGrey -g $i -s $(echo $i | sed -E 's#.*/(.*)\.fa#\1#') -o $(echo $i | sed -E 's#.*/(.*)\.fa#./\1/#') -t 50 -r agaricales -l /mnt/sdb/proj/Thibault/marasmius/db/repbase_all_2022-11.fasta -c yes -e yes \>$(echo $i | sed -E 's#.*/(.*)\.fa#\1.earlgrey.log#') 2\>\&1;done>earlgrey.genomes.sh
		nohup bash earlgrey.assemblies.sh >earlgrey.assemblies.log 2>&1 &
		nohup bash earlgrey.genomes.sh >earlgrey.genomes.sh 2>&1 &
#orthofinder
	conda create -n orthofinder orthofinder #pasted installation message in orthofinder.install
	mkdir phylo/orthofinder
	for i in $(ls ../../annotation/funannotate/*/ -d |grep -v test |grep -v db |grep -v old); do echo cp $i\annotate_results/*proteins.fa $(echo $i|cut -d'/' -f5|sed 's#$#.faa#');done|bash
	nohup orthofinder -f proteomes/ -t 48 -a 48 > orthofinder.log 2>orthofinder.err &
#ITS identification
	wget https://microbiology.se/sw/ITSx_1.1.3.tar.gz
	tar xvfz ITSx_1.1.3.tar.gz
	conda create -n itsx
	conda install hmmer #in itsx env
	mkdir its_id/
	nohup ../../ITSx_1.1.3/ITSx -i ../assemblies/6_cleaned/baemy1.fa -o baemy1 --cpu 20 --preserve T --save_regions ITS1,5_8S,ITS2 --partial 20 --nhmmer T >baemy1.itsx.log 2>&1 & #test
	for i in $(ls ../assemblies/6_cleaned/*.fa); do echo ../../ITSx_1.1.3/ITSx -i $i -o $(echo $i | sed -E 's#.*/([^/]+)\.fa#\1#') --cpu 20 --preserve T --save_regions ITS1,5_8S,ITS2 --partial 20 --nhmmer T \>$(echo $i | sed -E 's#.*/([^/]+)\.fa#\1.itsx.log#') 2\>\&1;done>itsx.sh
	nohup bash itsx.sh >itsx.log 2>&1 & #gymer1 doesn't have any full.fasta, trying with ITS1 instead
	for i in $(ls *full.fasta); do echo vsearch --derep_fulllength $i --output $(echo $i | sed 's#.fasta#.derep.fasta#') --sizeout --uc derep.uc \>$(echo $i | sed 's#.fasta#.derep.log#') 2\>\&1;done>vsearch.sh
	for i in $(ls *derep.fasta); do echo blastn -query $i -db ../db/ITS_eukaryote_sequences -outfmt \'6 qseqid sseqid pident length mismatch gapopen evalue bitscore\' -max_target_seqs 25 -perc_identity 85 -evalue 1e-20 -num_threads 50 \>$(echo $i | sed 's#.fasta#.its-euk.tsv#') 2\>$(echo $i | sed 's#.fasta#.its-euk.log#');done>blastn.sh #against blast ITS eukaryote
	makeblastdb -in sh_general_release_dynamic_19.02.2025.fasta -dbtype nucl #I also to it against the unite db (in ../db/)
	for i in $(ls *derep.fasta); do echo blastn -query $i -db ../db/sh_general_release_dynamic_19.02.2025.fasta -outfmt \'6 qseqid sseqid pident length mismatch gapopen evalue bitscore\' -max_target_seqs 25 -perc_identity 85 -evalue 1e-20 -num_threads 50 \>$(echo $i | sed 's#.fasta#.its-unite.tsv#') 2\>$(echo $i | sed 's#.fasta#.its-unite.log#');done>blastn-unite.sh
	for i in $(ls *.tsv); do sort -k3,3rn $i>$(echo $i|sed 's#.tsv#.sorted.tsv#');done
#gathering data
	#busco (~/Documents/WORK/MARASMIUS/busco/)
		mkdir assemblies
		scp thibault@emb-p-susrv01.emb.su.se:/mnt/sdb/proj/Thibault/marasmius/assemblies/5_funannotate/busco_results/BUSCO_*/short*json ./assemblies/
		mkdir genomes
		scp thibault@emb-p-susrv01.emb.su.se:/mnt/sdb/proj/Thibault/marasmius/genomes/busco_results/BUSCO_*/short*json .genomes/
		mkdir all_phylo
		scp thibault@emb-p-susrv01.emb.su.se:/mnt/sdb/proj/Thibault/marasmius/phylo/iqtree/busco/BUSCO_*/short*json .
		mkdir all-old
		scp thibault@emb-p-susrv01.emb.su.se:/mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/busco_results/BUSCO_*/short*json ./all-old/
		scp thibault@emb-p-susrv01.emb.su.se:/mnt/sdb/proj/Thibault/marasmius/assemblies/1_masurca/busco_results/BUSCO_*/short*json ./all-old/
		#imported on R, see busco.R
	#assembly features (~/Documents/WORK/MARASMIUS/assemblies/)
		scp thibault@emb-p-susrv01.emb.su.se:/mnt/sdb/proj/Thibault/marasmius/assemblies/6_cleaned/*fa .
		scp thibault@emb-p-susrv01.emb.su.se:/mnt/sdb/proj/Thibault/marasmius/genomes/cleaned/*fa .
		bash genome-assembly-stats.sh ./ > genome-assembly-stats.tsv
	#funannotate (~/Documents/WORK/MARASMIUS/annotation/funannotate)
		scp thibault@emb-p-susrv01.emb.su.se:/mnt/sdb/proj/Thibault/marasmius/annotation/funannotate/*/annotate_results/*gff3 . 
		bash genome-annotation-stats.sh ./ > genome-annotation-stats.tsv
#TO continue #
	#earlgrey (~/Documents/WORK/MARASMIUS/annotation/earlgrey)
		scp thibault@emb-p-susrv01.emb.su.se:/mnt/sdb/proj/Thibault/marasmius/annotation/earlgrey/*/*_EarlGrey/*_summaryFiles/*.familyLevelCount.txt .
		bash covmatrix.sh
		#import coverage_superfamily_matrix on R, see genome-stats.R
	#iqtree (~/Documents/WORK/MARASMIUS/phylo/iqtree )
		scp thibault@emb-p-susrv01.emb.su.se:/mnt/sdb/proj/Thibault/marasmius/phylo/iqtree/nt_90/supermatrix/SUPERMATRIX.phylip.cf.tree .
######
#END #
######

#other decontamination tools to think about: https://github.com/h836472/ContScout https://github.com/ncbi/fcs
#scaffolding tools: ragtag, yaHS

#### OTHER COMMANDS NOT SELECTED FOR FURTHER STEPS ####
#decontaminate reads
	conda create -n kraken2 kraken2
	cd data
	mkdir kraken
	cd kraken
	kraken2-build --standard --db kraken2_db --threads 64 >build.log 2>&1 & #completed database >300Go, trying a lighter one
	wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20250402.tar.gz	& #standard 90Gb
	wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20250402.tar.gz & #Standard DB capped at 16 GB
	wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20250402.tar.gz & #Standard plus Refeq protozoa, fungi & plant, capped at 16 GB
	kraken2-build --download-taxonomy --db k2_pluspf16
	kraken2-build --download-taxonomy --db k2_standard16
	tar -xvzf k2_pluspfp_16gb_20250402.tar.gz -C k2_pluspf16
	tar -xvzf k2_standard_16gb_20250402.tar.gz -C k2_standard16
	#kraken need to have a taxid in the fasta heaer, remplacing them
	cd ../../genomes/
	sed 's/^>.*/>137718/' colco1.raw.sort.masked.fa > colco1.raw.sort.masked_taxid.fa
	sed 's/^>.*/>38945/' flave1.raw.sort.masked.fa > flave1.raw.sort.masked_taxid.fa 
	sed 's/^>.*/>762847/' gyman1.raw.sort.masked.fa > gyman1.raw.sort.masked_taxid.fa 
	sed 's/^>.*/>564482/' gymlu1.raw.sort.masked.fa > gymlu1.raw.sort.masked_taxid.fa 
	sed 's/^>.*/>5353/' lened1.raw.sort.masked.fa > lened1.raw.sort.masked_taxid.fa 
	sed 's/^>.*/>585013/' marcr1.raw.sort.masked.fa > marcr1.raw.sort.masked_taxid.fa 
	sed 's/^>.*/>1328755/' marfi1.raw.sort.masked.fa > marfi1.raw.sort.masked_taxid.fa 
	sed 's/^>.*/>181124/' maror1.raw.sort.masked.fa > maror1.raw.sort.masked_taxid.fa 
	sed 's/^>.*/>3065222/' marsp1.raw.sort.masked.fa > marsp1.raw.sort.masked_taxid.fa 
	sed 's/^>.*/>2682957/' marsc1.raw.sort.masked.fa > marsc1.raw.sort.masked_taxid.fa 
	sed 's/^>.*/>585030/' marte1.raw.sort.masked.fa > marte1.raw.sort.masked_taxid.fa 
	sed 's/^>.*/>221102/' monro1.raw.sort.masked.fa > monro1.raw.sort.masked_taxid.fa 
	sed 's/^>.*/>5334/' schco1.raw.sort.masked.fa > schco1.raw.sort.masked_taxid.fa
	cd ../data/kraken/ #add them to db (same for others + in k2_standard16)
	kraken2-build --add-to-library ../../genomes/colco1.raw.sort.masked_taxid.fa --db k2_pluspf16 #k2_standard16
	kraken2-build --build --db k2_standard16 --threads 32 >k2s16_build.log 2>&1 &
	kraken2-build --build --db k2_pluspf16 --threads 32 >k2pf16_build.log 2>&1 &
	kraken2 --db k2_standard16 --threads 16 --report baemy1_s16_report.txt --classified-out baemy1_s16_classified-out.fq --unclassified-out baemy1_s16_unclassified-out.fq --gzip-compressed --output baemy1_s16_clean.txt ../ccsreads/baemy1_004/pr_017_004.ccsreads.fastq.gz >baemy1_k2s16.log 2>&1 & #same for others)
	#not very better, I give up and go back to haplotype 1 from hifiasm and blobtools
#blob
		#generate hit files
             wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz #BLAST nt.gz db is 378G (compressed), online instead
             nohup blastn -query assemblies/hifiasm/pr_017_002.bp.p_ctg.gfa.fa -db nt -remote -outfmt 6 >hit/pr_017_002.tsv 2>hit/pr_017_002.err & [1] 3805612 #same for others, too long
             #2nd trial splitting the genomes
                 seqkit split -p 96 ../../assemblies/masurca/gymer/CA/primary.genome.scf.fasta #splitting command for the two illumina genomes
                 cd hit
                 mkdir pr_017_002 #same for others
                 vim pr_017002.sh
                     #!/bin/bash
                     eval "$(conda shell.bash hook)"
                     for i in $(grep ">" /mnt/sda/proj/Thibault/marasmius/assemblies/hifiasm/pr_017_002.bp.p_ctg.gfa.fa | cut -f1 -d" " | sed 's/>//'); do
                             conda activate samtools;
                             samtools faidx /mnt/sda/proj/Thibault/marasmius/assemblies/hifiasm/pr_017_002.bp.p_ctg.gfa.fa $i -o /mnt/sda/proj/Thibault/marasmius/hit/pr_017_002/$i\.fa;
                             conda deactivate
                                 conda activate blast
                                 blastn -query /mnt/sda/proj/Thibault/marasmius/hit/pr_017_002/$i\.fa -db nt -outfmt 6 -remote -out /mnt/sda/proj/Thibault/marasmius/hit/pr_017_002/$i\_nt.tsv;
                             conda deactivate;
                             rm /mnt/sda/proj/Thibault/marasmius/hit/pr_017_002/$i\.fa
                         done
                 sed 's/002/003/g' pr_017_002.sh > pr_017_003.sh #same for others
                 chmod 755 *.sh
                 nohup ./pr_017_002.sh > 002.logerr 2>&1 & [1] 3835926 #ok for 8 and 11, too long and killed for the others
                 #concatenation (only for 008 and 011)
                     cd hit/pr_017_008/
                     for i in $(ls); do cat $i;done > pr_017_008-nt.tsv
                     rm ptg*
             #3rd trial with the nt database in spore (same for others)
                 nohup blastn -query assemblies/hifiasm/pr_017_002.bp.p_ctg.gfa.fa -db /mnt/sda/blast_databases/nt -outfmt 6 -out hit/pr_017_002/pr_017_002-nt.tsv -num_threads 9 >hit/pr_017_002/pr_017_002.logerr 2>&1 & [1] 4011000 #too long
#assembly
	#canu
		conda create -n canu canu
		mkdir assemblies/canu
		for i in $(ls ccsreads/ | grep pr); do echo canu -p $i -d assemblies/canu/$i genomeSize=60m -pacbio-hifi ccsreads/$i/$i\.ccsreads.fastq.gz \>assemblies/canu/$i.stdouterr 2\>\&1;done | bash &
		#005 and 009 show low coverage (8.58 and 3.62), rescue with 48m and 22m genome size, 009 and 011 are doomed (0.72 and 0.61) 
		canu -p pr_017_005 -d assemblies/canu/pr_017_005/run48m genomeSize=48m -pacbio-hifi ccsreads/pr_017_005/pr_017_005.ccsreads.fastq.gz >assemblies/canu/pr_017_005/pr_017_005-run48m.stdouterr 2>&1 &
		canu -p pr_017_009 -d assemblies/canu/pr_017_009/run20m genomeSize=20m -pacbio-hifi ccsreads/pr_017_009/pr_017_009.ccsreads.fastq.gz >assemblies/canu/pr_017_009/pr_017_009-run20m.stdouterr 2>&1 &
		#busco
			cd assemblies/canu
			mkdir busco_results
			busco -i pr_017_002/pr_017_002.contigs.fasta -m geno -l basidiomycota_odb10 -c 12 --out_path busco_results -q >busco_results/busco02.outerr 2>&1 & #same for others
			busco -i pr_017_005/run48m/pr_017_005.contigs.fasta -m geno -l basidiomycota_odb10 -c 12 --out_path busco_results -q >busco_results/busco05.outerr 2>&1 &
			busco -i pr_017_009/run20m/pr_017_009.contigs.fasta -m geno -l basidiomycota_odb10 -c 12 --out_path busco_results -q >busco_results/busco09.outerr 2>&1 &
			for i in $(ls | grep BUSCO); do sed -n 9p $i/*.txt; done > busco_results.txt #formatting in calc
		#quast
			cd assemblies/canu
			mkdir quast_results/
			quast pr_017_002/pr_017_002.contigs.fasta -o quast_results/pr_017_002 >quast_results/quast02.outerr 2>&1 & #same for others
			quast pr_017_005/run48m/pr_017_005.contigs.fasta -o quast_results/pr_017_005 >quast_results/quast05.outerr 2>&1 &
			quast pr_017_009/run20m/pr_017_009.contigs.fasta -o quast_results/pr_017_009 >quast_results/quast09.outerr 2>&1 &
			cd quast_results/
			for i in $(ls | grep pr); do tail -n 1 $i/transposed_report.tsv;done > transposed_report.tsv #paste header and transpose via calc
	#HiFlye
		conda create -n flye flye
		mkdir assemblies/flye
		for i in $(ls ccsreads | grep pr); do flye --pacbio-hifi ccsreads/$i/$i\.ccsreads.fastq.gz -t 10 -o assemblies/flye/$i; done >assemblies/flye/flye.outerr 2>&1 &
		for i in $(ls| grep pr); do mv $i/assembly.fasta $i/$i\_flye-assembly.fasta;done #rename
		#busco
			cd assemblies/flye
			mkdir busco_results
			nohup bash -c 'for i in $(ls | grep pr); do busco -i $i/$i\_flye-assembly.fasta -m geno -l basidiomycota_odb10 -c 88 --out_path busco_results -q; done' >busco_results/nohup.outerr 2>&1 &
			#08 and 011 faild, redo without quiet mode to see
			busco -i pr_017_008/pr_017_008_flye-assembly.fasta -m geno -l basidiomycota_odb10 -c 88 --out_path busco_results > busco_results/busco08.out 2>&1 & [1] 2171436 #no match found
			busco -i pr_017_011/pr_017_011_flye-assembly.fasta -m geno -l basidiomycota_odb10 -c 88 --out_path busco_results > busco_results/busco11.out 2>&1 &[1] 2211693 #no match found
			cd busco_results/
			for i in $(ls | grep BUSCO); do sed -n 9p $i/*.txt; done > busco_results.txt #formatting in calc
		#quast
			cd assemblies/flye
			mkdir quast_results/
			for i in $(ls | grep pr); do echo quast $i/$i\_flye-assembly.fasta -o quast_results/$i; done | bash >quast_results/quast.out &
			cd quast_results/
			for i in $(ls | grep pr); do tail -n 1 $i/transposed_report.tsv;done > transposed_report.tsv #paste header and transpose via calc
	#Peregrine
		conda create -n peregrine bioconda::peregrine-2021
		conda activate peregrine
		for i in $(ls ccsreads | grep pr); do pg_asm ccsreads/$i/$i\.ccsreads.fastq.gz assemblies/peregrine/$i 48 1 ;done >>assemblies/peregrine/peregrine.outerr 2>&1 &
		#error: thread 'main' panicked at 'called `Result::unwrap()` on an `Err` value: Error { kind: InvalidData, message: "stream did not contain valid UTF-8" }', src/bin/utils/build_sdb.rs:249:35
		#try with one thread
		for i in $(ls ccsreads | grep pr); do pg_asm ccsreads/$i/$i\.ccsreads.fastq.gz assemblies/peregrine/$i 1 1 ;done >>assemblies/peregrine/peregrine.outerr 2>&1 &
		#same, try with absolute path
		for i in $(ls ccsreads | grep pr); do pg_asm /mnt/sda/proj/Thibault/raw-data/pr_017/ccsreads/$i/$i\.ccsreads.fastq.gz /mnt/sda/proj/Thibault/raw-data/pr_017/assemblies/peregrine/$i 1 1 ;done >>assemblies/peregrine/peregrine-abs.outerr 2>&1 &
		#same, try with a text file containning the sequence paths, as showed in the tutorial https://github.com/cschin/peregrine-2021
		for i in $(ls ccsreads | grep pr); do find /mnt/sda/proj/Thibault/raw-data/pr_017/ccsreads/$i/ | grep gz > assemblies/peregrine/$i\-reads.lst;done
		for i in $(ls ccsreads | grep pr); do pg_asm /mnt/sda/proj/Thibault/raw-data/pr_017/assemblies/peregrine/$i\-reads.lst /mnt/sda/proj/Thibault/raw-data/pr_017/assemblies/peregrine/$i 1 1 ;done >>assemblies/peregrine/peregrine-lst.outerr 2>&1 &
		#success! redo with 48 threads
		kill 2295345
		for i in $(ls ccsreads | grep pr); do pg_asm /mnt/sda/proj/Thibault/raw-data/pr_017/assemblies/peregrine/$i\-reads.lst /mnt/sda/proj/Thibault/raw-data/pr_017/assemblies/peregrine/$i 48 1 ;done >>assemblies/peregrine/peregrine-final.outerr 2>&1 &
		#its only taking one CPU, try with 2 "partitions" to see
		kill 2295477
		for i in $(ls ccsreads | grep pr); do pg_asm /mnt/sda/proj/Thibault/raw-data/pr_017/assemblies/peregrine/$i\-reads.lst /mnt/sda/proj/Thibault/raw-data/pr_017/assemblies/peregrine/$i 4 4 ;done >>assemblies/peregrine/peregrine-4t4p.outerr 2>&1 &
		#taking 4 CPU, now launch it with more
		kill 2295885
		for i in $(ls ccsreads | grep pr); do pg_asm /mnt/sda/proj/Thibault/raw-data/pr_017/assemblies/peregrine/$i\-reads.lst /mnt/sda/proj/Thibault/raw-data/pr_017/assemblies/peregrine/$i 2 48 ;done >>assemblies/peregrine/peregrine-final-final.outerr 2>&1 &
		#now taking 2 or more CPU, try with 8 threads and 8 partitions
		kill 2295999
		for i in $(ls ccsreads | grep pr); do pg_asm /mnt/sda/proj/Thibault/raw-data/pr_017/assemblies/peregrine/$i\-reads.lst /mnt/sda/proj/Thibault/raw-data/pr_017/assemblies/peregrine/$i 8 8 ;done >>assemblies/peregrine/peregrine-final-final-final.outerr 2>&1 &
		grep "pg_asm run end" assemblies/peregrine/peregrine-final-final-final.outerr #check everyone is done
		for i in $(ls | grep pr | grep -v reads); do mv $i/asm_ctgs_m_p.fa $i/$i\_peregrine-assembly.fa;done #rename
		#busco
			cd assemblies/peregrine
			mkdir assemblies/peregrine/busco_results
			nohup bash -c 'for i in $(ls | grep pr | grep -v reads); do busco -i $i/$i\_peregrine-assembly.fa -m geno -l basidiomycota_odb10 -c 88 --out_path busco_results -q; done' >busco_results/nohup.outerr 2>&1 &
			cd busco_results/
			for i in $(ls | grep BUSCO); do sed -n 9p $i/*.txt; done > busco_results.txt #formatting in calc
		#quast
			cd assemblies/peregrine
			mkdir quast_results/
			for i in $(ls | grep pr | grep -v reads); do quast $i/$i\_peregrine-assembly.fa -o quast_results/$i; done >quast_results/quast.out &
			cd quast_results/
			for i in $(ls | grep pr); do tail -n 1 $i/transposed_report.tsv;done > transposed_report.tsv #paste header and transpose via calc
	#nextDenovo
		conda create -n nextdenovo nextdenovo
		conda activate nextdenovo
		mkdir assemblies/nextdenovo
		#make the cfg files
		#copy, annotate and edit (default unless genome size = 60m) a model file in assemblies/nextdenovo/run.cfg from https://nextdenovo.readthedocs.io/en/latest/OPTION.html
		for i in $(ls ccsreads | grep pr); do sed -e "s|input\.fofn|/mnt/sda/proj/Thibault/raw-data/pr_017/assemblies/peregrine/$i-reads.lst|" assemblies/nextdenovo/run.cfg | sed "s|workdir/| /mnt/sda/proj/Thibault/raw-data/pr_017/assemblies/nextdenovo/$i|" > assemblies/nextdenovo/$i.cfg; done #change input (using the same file as peregrine) and workdir
		cd assemblies/nextdenovo
		for i in $(ls | grep pr); do nextDenovo $i;done >nextdenovo.outerr 2>&1 &
		grep -A2 Suggested nextdenovo.outerr #5, 8, 9 and 11 are not complete, not enough coverage, let's do another round with 40m
		cp pr_017_005.cfg pr_017_005-40m.cfg #for those four, copy and manually change genome size from 60m to 40m
		rm -r pr_017_005/ #for those four, remove the folder (otherwise it uses the old files and old genome size)
		nextDenovo pr_017_005-40m.cfg >nextdenovo05-40m.outerr 2>&1 & #ok
		nextDenovo pr_017_008-40m.cfg >nextdenovo08-40m.outerr 2>&1 & #not working again, trying with 20m
		nextDenovo pr_017_008-20m.cfg >nextdenovo08-20m.outerr 2>&1 & #not working, I give up
		nextDenovo pr_017_009-40m.cfg >nextdenovo09-40m.outerr 2>&1 & #not working again, trying with 20m
		nextDenovo pr_017_009-20m.cfg >nextdenovo09-20m.outerr 2>&1 & #ok
		nextDenovo pr_017_011-40m.cfg >nextdenovo11-40m.outerr 2>&1 & #not working again, trying with 20m
		nextDenovo pr_017_011-20m.cfg >nextdenovo11-20m.outerr 2>&1 & #not working, I give up
		for i in $(ls | grep pr | grep -v cfg); do mv $i/03.ctg_graph/nd.asm.fasta $i/03.ctg_graph/$i\_nextdenovo-assembly.fa;done #rename
		#busco
			mkdir busco_results
			nohup bash -c 'for i in $(ls | grep pr | grep -v cfg); do busco -i $i/03.ctg_graph/$i\_nextdenovo-assembly.fa -m geno -l basidiomycota_odb10 -c 88 --out_path busco_results -q; done' >busco_results/nohup.outerr 2>&1 &
			for i in $(ls | grep BUSCO); do sed -n 9p $i/*.txt; done > busco_results.txt #formatting in calc
		#quast
			cd assemblies/nextdenovo
			mkdir quast_results/
			for i in $(ls | grep pr | grep -v cfg); do quast $i/03.ctg_graph/$i\_nextdenovo-assembly.fa -o quast_results/$i; done >quast_results/quast.out &
			cd quast_results/
			for i in $(ls | grep pr); do tail -n 1 $i/transposed_report.tsv;done > transposed_report.tsv #paste header and transpose via calc
	#falcon
		conda create -n falcon pb-assembly
		mkdir assemblies/falcon
		#looks complicated to run, cancelled
		#Spades (M1 Gymer 001 and M21 Marro 012)
			conda create -n spades spades
			conda activate spades
			mkdir assemblies/spades
			mkdir assemblies/spades/gymer
			spades.py -1 /mnt/sda/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M1-Gymer_S112-S113_L003_R1_001.fastq.gz -2 /mnt/sda/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M1-Gymer_S112-S113_L003_R2_001.fastq.gz -t 88 -m 502 -o /mnt/sda/proj/Thibault/marasmius/assemblies/spades/gymer/ > /mnt/sda/proj/Thibault/marasmius/assemblies/spades/gymer/spades.logerr 2>&1 &
			mkdir assemblies/spades/marro
			nohup spades.py -1 /mnt/sda/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M21-Marro_S114-S115_L003_R1_001.fastq.gz -2 /mnt/sda/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M21-Marro_S114-S115_L003_R2_001.fastq.gz -t 16 -m 200 -o /mnt/sda/proj/Thibault/marasmius/assemblies/spades/marro/ > /mnt/sda/proj/Thibault/marasmius/assemblies/spades/marro/spades.logerr 2>&1 &
			#busco (same for both)
				cd spades/gymer
				mkdir busco_results
				busco -i scaffolds.fasta -m geno -l basidiomycota_odb10 -c 88 --out_path busco_results -q >busco_results/busco.logerr 2>&1 & [1] 3971503 #killed, RAM consumption exceeded with 88 and 24 CPU, try again with 12 CPUs
				busco -i scaffolds.fasta -m geno -l basidiomycota_odb10 -c 12 --out_path busco_results -q >busco_results/busco.logerr 2>&1 & [10] 4011740
				busco -i scaffolds.fasta -m geno -l basidiomycota_odb10 -c 24 --out_path busco_results -q >busco_results/busco.logerr 2>&1 &
				sed -n 9p busco_results/BUSCO_scaffolds.fasta/short_summary.specific.basidiomycota_odb10.BUSCO_scaffolds.fasta.txt #formatting in calc
			#quast (same for both)
				cd spades/gymer/
				mkdir quast_results/
				quast scaffolds.fasta -o quast_results/ >quast_results/quast.out & 
				quast scaffolds.fasta -o quast_results/ >quast_results/quast.out &
				tail -n 1 quast/results/transposed_report.tsv #pasted locally with header
	#a lot of duplications found in 002, 003, 004, 006, and 007, investigation
		mkdir assemblies/duplications
		cp -r busco_results/BUSCO_mycal1.clean.sort.masked.fa/run_agaricales_odb10/busco_sequences/multi_copy_busco_sequences/*fna duplications/mycal1_multi_busco
		cp -r busco_results/BUSCO_rhodo1.clean.sort.masked.fa/run_agaricales_odb10/busco_sequences/multi_copy_busco_sequences/*fna duplications/rhodo1_multi_busco
		cp -r busco_results/BUSCO_baemy1.clean.sort.masked.fa/run_agaricales_odb10/busco_sequences/multi_copy_busco_sequences/*fna duplications/baemy1_multi_busco
		cp -r busco_results/BUSCO_mycsc1.clean.sort.masked.fa/run_agaricales_odb10/busco_sequences/multi_copy_busco_sequences/*fna duplications/mycsc1_multi_busco 
		cp -r busco_results/BUSCO_mycal2.clean.sort.masked.fa/run_agaricales_odb10/busco_sequences/multi_copy_busco_sequences/*fna duplications/mycal2_multi_busco
		#bring haplotypes (same for others)
		cp hifiasm/pr_017_002.bp.hap1.p_ctg.gfa duplications/mycal1.hap1.gfa
		cp hifiasm/pr_017_002.bp.hap2.p_ctg.gfa duplications/mycal1.hap2.gfa
	#gfa to fasta (same for others)
		awk '/^S/{print ">"$2"\n"$3}' baemy1.hap1.gfa | fold > baemy1.hap1.fa	
	#homology alignment
		#mummer
			cd duplications
			mummer -maxmatch -l 1000 ../funannotate/run1/baemy1.clean.sort.masked.fa ../funannotate/run1/baemy1.clean.sort.masked.fa >baemy1-l1000max.delta
			mummerplot -R ../funannotate/run1/baemy1.clean.sort.masked.fa -Q ../funannotate/run1/baemy1.clean.sort.masked.fa -p baemy1-l1000max --png baemy1-l1000max.delta
	#look at family classification
		cd funannotate/run1/busco_results/BUSCO_baemy1.clean.sort.masked.fa/run_agaricales_odb10
		grep "Duplicated" full_table.tsv | cut -f 10 | sort | uniq -c | sort -nr | sed -E 's/^\s+//; s/\s/\t/' | head #same with complete genes, compare
	#look at distance to scaffold edge
		grep "Complete" full_table.tsv | cut -f 3,4,5 | awk -v OFS='\t' '{print $1, sqrt(($2-$3)*($2-$3))}' > rhodo1-complete.pos #get central position of the gene
		bioawk -c fastx '{print $name length($seq)}' ../../../rhodo1.clean.sort.masked.fa #get scaffold length
		awk 'NR==FNR {scaffold_length[$1] = $2;next}{if ($1 in scaffold_length) {left = $2;right = scaffold_length[$1] - $2;min_distance = (left < right) ? left : right;print min_distance;}}' rhodo-scaff-length rhodo1-complete.pos > rhodo-complete-distance-edge #same for duplicated
	#look at the duplication numbers
		cd funannotate/run1/busco_results/BUSCO_mycal2.clean.sort.masked.fa/run_agaricales_odb10
		grep "Duplicated" full_table.tsv | cut -f1 | uniq -c | sort | cut -f7 -d' ' | uniq -c
	#busco against haplotypes
		cd duplications
		busco -i baemy1.hap1.fa -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 50 --out_path busco_results -q >busco_results/baemy1.hap1.busco.log 2>&1 & 
		busco -i baemy1.hap2.fa -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 50 --out_path busco_results -q >busco_results/baemy1.hap2.busco.log 2>&1 & 
		busco -i mycal1.hap1.fa -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 50 --out_path busco_results -q >busco_results/mycal1.hap1.busco.log 2>&1 & 
		busco -i mycal1.hap2.fa -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 50 --out_path busco_results -q >busco_results/mycal1.hap2.busco.log 2>&1 & 
		busco -i mycal2.hap1.fa -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 50 --out_path busco_results -q >busco_results/mycal2.hap1.busco.log 2>&1 & 
		busco -i mycal2.hap2.fa -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 50 --out_path busco_results -q >busco_results/mycal2.hap2.busco.log 2>&1 & 
		busco -i mycsc1.hap1.fa -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 50 --out_path busco_results -q >busco_results/mycsc1.hap1.busco.log 2>&1 & 
		busco -i mycsc1.hap2.fa -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 50 --out_path busco_results -q >busco_results/mycsc1.hap2.busco.log 2>&1 & 
		busco -i rhodo1.hap1.fa -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 50 --out_path busco_results -q >busco_results/rhodo1.hap1.busco.log 2>&1 & 
		busco -i rhodo1.hap2.fa -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 50 --out_path busco_results -q >busco_results/rhodo1.hap2.busco.log 2>&1 & 
		#ghostx
			conda create -n ghostx -b bioconda ghostx
			#database (same for others)
				ghostx db -i baemy1.hap1.fa -o baemy1.hap1.db -t d > baemy1.hap1.db.log 2>&1 &
				ghostx db -i baemy1.hap2.fa -o baemy1.hap2.db -t d > baemy1.hap2.db.log 2>&1 &
			#align
				for i in  $(ls baemy1_multi_busco/*); do ghostx aln -i $i -o $(echo $i | sed 's:fna:hap1\.aln:') -d baemy1.hap1.db -q d -a 10;done > baemy1.hap1.aln.log 2>&1 & 
				for i in  $(ls baemy1_multi_busco/*); do ghostx aln -i $i -o $(echo $i | sed 's:fna:hap2\.aln:') -d baemy1.hap2.db -q d -a 10;done > baemy1.hap2.aln.log 2>&1 & 
				for i in  $(ls mycal1_multi_busco/*); do ghostx aln -i $i -o $(echo $i | sed 's:fna:hap1\.aln:') -d mycal1.hap1.db -q d -a 10;done > mycal1.hap1.aln.log 2>&1 & 
				for i in  $(ls mycal1_multi_busco/*); do ghostx aln -i $i -o $(echo $i | sed 's:fna:hap2\.aln:') -d mycal1.hap2.db -q d -a 10;done > mycal1.hap2.aln.log 2>&1 & 
				for i in  $(ls mycal2_multi_busco/*); do ghostx aln -i $i -o $(echo $i | sed 's:fna:hap1\.aln:') -d mycal2.hap1.db -q d -a 10;done > mycal2.hap1.aln.log 2>&1 & 
				for i in  $(ls mycal2_multi_busco/*); do ghostx aln -i $i -o $(echo $i | sed 's:fna:hap2\.aln:') -d mycal2.hap2.db -q d -a 10;done > mycal2.hap2.aln.log 2>&1 & 
				for i in  $(ls mycsc1_multi_busco/*); do ghostx aln -i $i -o $(echo $i | sed 's:fna:hap1\.aln:') -d mycsc1.hap1.db -q d -a 10;done > mycsc1.hap1.aln.log 2>&1 & 
				for i in  $(ls mycsc1_multi_busco/*); do ghostx aln -i $i -o $(echo $i | sed 's:fna:hap2\.aln:') -d mycsc1.hap2.db -q d -a 10;done > mycsc1.hap2.aln.log 2>&1 & 
				for i in  $(ls rhodo1_multi_busco/*); do ghostx aln -i $i -o $(echo $i | sed 's:fna:hap1\.aln:') -d rhodo1.hap1.db -q d -a 10;done > rhodo1.hap1.aln.log 2>&1 & 
				for i in  $(ls rhodo1_multi_busco/*); do ghostx aln -i $i -o $(echo $i | sed 's:fna:hap2\.aln:') -d rhodo1.hap2.db -q d -a 10;done > rhodo1.hap2.aln.log 2>&1 & 
			#count how many in hap1, 2, both or none
			for i in *fna; do h1="${i/.fna/.hap1.aln}";h2="${i/.fna/.hap2.aln}";wc1=$(wc -l < "$h1");wc2=$(wc -l < "$h2");if [ "$wc1" -eq 0 ] && [ "$wc2" -eq 0 ] ; then echo "0";elif [ "$wc1" -ne 0 ] && [ "$wc2" -eq 0 ] ; then echo "1";elif [ "$wc1" -eq 0 ] && [ "$wc2" -ne 0 ] ; then echo "2";else echo "3";fi; done | sort | uniq -c
	#checking common duplicated genes between mycal1, mycal2 and mycsc1
		cd assemblies/funannotate/run1/busco_results
		grep "Duplicated" BUSCO_mycal1.clean.sort.masked.fa/run_agaricales_odb10/full_table.tsv | cut -f 1,2,9,10 | uniq > mycal1.duplicated #wc -l=348
		grep "Duplicated" BUSCO_mycal2.clean.sort.masked.fa/run_agaricales_odb10/full_table.tsv | cut -f 1,2,9,10 | uniq > mycal2.duplicated #wc -l=920
		grep "Duplicated" BUSCO_mycsc1.clean.sort.masked.fa/run_agaricales_odb10/full_table.tsv | cut -f 1,2,9,10 | uniq > mycsc1.duplicated #wc -l=590
		grep -Ff mycal2.duplicated mycal1.duplicated > mycal1-2.duplicated #wc -l=171
		grep -Ff mycal2.duplicated mycsc1.duplicated > mycal2-mycsc.duplicated #wc -l=151
		grep -Ff mycal1.duplicated mycsc1.duplicated > mycal1-mycsc.duplicated #wc -l=60
		grep -Ff mycal1-2.duplicated mycal2-mycsc.duplicated > mycal1-2-mycsc.duplicated #wc -l=41 >all the rest follows
#purge_dups
#trying another round for baemy
		cd 3_purge_dups/baemy1
		mkdir run2
		cp purged.fa baemy1.hifiasm.blob.purged.fa #purged.fa has to be the output name
		pd_config.py -n baemy1.config.json -l /mnt/sdb/proj/Thibault/marasmius/assemblies/3_purge_dups/baemy1/run2 /mnt/sdb/proj/Thibault/marasmius/assemblies/3_purge_dups/baemy1/baemy1.hifiasm.blob.purged.fa /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/baemy1_004/pr_017_004.ccsreads.fastq.gz
		vim baemy1.config.json #add lineage agaricales
		run_purge_dups.py /mnt/sdb/proj/Thibault/marasmius/assemblies/3_purge_dups/baemy1/baemy1.config.json /mnt/sdb/proj/Thibault/miniconda3/envs/purge_dups/bin baemy1 > baemy1.purge_dups.log 2>&1 & #still trying the wrapped command, but same error, switch to full pipeline script again
		cd run2/
		cp ../baemy1_purge-dups.sh .
		sed -i 's/blob.fa/blob.purged.fa/g' baemy1_purge-dups.sh
		nohup bash baemy1_purge-dups.sh 2>&1 & #error at pbcstat command, cutoffs at 5,0,0,0,0, PB.cov.wig empty, it's probably hopeless but trying with a lower min secondary-to-primary score ratio alignment (from one github issue), here are the command different from initial
		minimap2 -xasm20 /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/baemy1_004/pr_017_004.ccsreads.fastq.gz baemy1.hifiasm.blob.fa -t 40 -p 0.2 | gzip -c - > baemy1.paf.gz #-p0.2
		minimap2 -xasm5 -DP baemy1.hifiasm.blob.purged.split baemy1.hifiasm.blob.purged.split -t 40 -p 0.2 | gzip -c - > baemy1.hifiasm.blob.purged.fa.split.self.paf.gz #-p 0.2
		get_seqs dups.bed baemy1.hifiasm.blob.fa #no -e option
		#removed 4 scaffolds/90kb, but no difference in gene content, I am trying on mycal2 which has more duplicated genes to see if it can make a difference, same procedure
		#all is removed with a p 0.2 instead of 0.8, same with 0.5, 0.7, same with 0.75
	#purge_haplotigs
		cd assemblies/
		mkdir purge_haplotigs
		nohup minimap2 -a -t 5 -x map-hifi /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/baemy1.bp.p_ctg.gfa.fa /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/baemy1_004/pr_017_004.ccsreads.fastq.gz -o baemy1_hifi.sam > baemy1_hifi.minimap2 2>&1 & 
		nohup minimap2 -a -t 5 -x map-hifi /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/colpe1.bp.p_ctg.gfa.filtered.fna /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/colpe1_005/pr_017_005.ccsreads.fastq.gz -o colpe1_hifi.sam > colpe1_hifi.minimap2 2>&1 &
		nohup minimap2 -a -t 5 -x map-hifi /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/gymaq1.bp.p_ctg.gfa.filtered.fna /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/gymaq1_010/pr_017_010.ccsreads.fastq.gz -o gymaq1_hifi.sam > gymaq1_hifi.minimap2 2>&1 & 
		nohup minimap2 -a -t 5 -x map-hifi /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/marsi1.bp.p_ctg.gfa.fa /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/marsi1_009/pr_017_009.ccsreads.fastq.gz -o marsi1_hifi.sam > marsi1_hifi.minimap2 2>&1 &
		nohup minimap2 -a -t 5 -x map-hifi /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/mycal1.bp.p_ctg.gfa.filtered.fna /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/mycal1_002/pr_017_002.ccsreads.fastq.gz -o mycal1_hifi.sam > mycal1_hifi.minimap2 2>&1 & 
		nohup minimap2 -a -t 5 -x map-hifi /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/mycal2.bp.p_ctg.gfa.fa /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/mycal2_007/pr_017_007.ccsreads.fastq.gz -o mycal2_hifi.sam > mycal2_hifi.minimap2 2>&1 & 
		nohup minimap2 -a -t 5 -x map-hifi /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/mycsc1.bp.p_ctg.gfa.filtered.fna /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/mycsc1_006/pr_017_006.ccsreads.fastq.gz -o mycsc1_hifi.sam > mycsc1_hifi.minimap2 2>&1 &
		nohup minimap2 -a -t 5 -x map-hifi /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/rhodo1.bp.p_ctg.gfa.filtered.fna /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/rhodo1_003/pr_017_003.ccsreads.fastq.gz -o rhodo1_hifi.sam > rhodo1_hifi.minimap2 2>&1 & 
		#trying with classic PacBio read option
		nohup minimap2 -a -t 5 -x map-pb /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/baemy1.bp.p_ctg.gfa.fa /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/baemy1_004/pr_017_004.ccsreads.fastq.gz -o baemy1_pb.sam > baemy1_pb.minimap2 2>&1 &
		nohup minimap2 -a -t 5 -x map-pb /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/colpe1.bp.p_ctg.gfa.filtered.fna /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/colpe1_005/pr_017_005.ccsreads.fastq.gz -o colpe1_pb.sam > colpe1_pb.minimap2 2>&1 &
		nohup minimap2 -a -t 5 -x map-pb /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/gymaq1.bp.p_ctg.gfa.filtered.fna /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/gymaq1_010/pr_017_010.ccsreads.fastq.gz -o gymaq1_pb.sam > gymaq1_pb.minimap2 2>&1 &
		nohup minimap2 -a -t 5 -x map-pb /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/marsi1.bp.p_ctg.gfa.fa /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/marsi1_009/pr_017_009.ccsreads.fastq.gz -o marsi1_pb.sam > marsi1_pb.minimap2 2>&1 &
		nohup minimap2 -a -t 5 -x map-pb /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/mycal1.bp.p_ctg.gfa.filtered.fna /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/mycal1_002/pr_017_002.ccsreads.fastq.gz -o mycal1_pb.sam > mycal1_pb.minimap2 2>&1 & 
		nohup minimap2 -a -t 5 -x map-pb /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/mycal2.bp.p_ctg.gfa.fa /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/mycal2_007/pr_017_007.ccsreads.fastq.gz -o mycal2_pb.sam > mycal2_pb.minimap2 2>&1 & 
		nohup minimap2 -a -t 5 -x map-pb /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/mycsc1.bp.p_ctg.gfa.filtered.fna /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/mycsc1_006/pr_017_006.ccsreads.fastq.gz -o mycsc1_pb.sam > mycsc1_pb.minimap2 2>&1 &
		nohup minimap2 -a -t 5 -x map-pb /mnt/sdb/proj/Thibault/marasmius/assemblies/1_hifiasm/rhodo1.bp.p_ctg.gfa.filtered.fna /mnt/sdb/proj/Thibault/marasmius/data/ccsreads/rhodo1_003/pr_017_003.ccsreads.fastq.gz -o rhodo1_pb.sam > rhodo1_pb.minimap2 2>&1 &
		#illumina
		cd assemblies/1_masurca/gymer/CA/
		bowtie2-build primary.genome.scf.filtered.fna gymer.filtered >gymer.bowtie2-build.log 2>&1 &
		cd ../../marro/CA/
		bowtie2-build primary.genome.scf.filtered.fna marro.filtered >marro.bowtie2-build.log 2>&1 &
		cd assemblies/3.2_purge_haplotigs/mapping/
		nohup bowtie2-align-s -x /mnt/sdb/proj/Thibault/marasmius/assemblies/1_masurca/gymer/CA/gymer.filtered -1 /mnt/sdb/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M1-Gymer_S112-S113_L003_R1_001.fastq.gz -2 /mnt/sdb/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M1-Gymer_S112-S113_L003_R2_001.fastq.gz -S gymer1.sam -p 16 > gymer1-bowtie2.log 2>&1 &
		nohup bowtie2-align-s -x /mnt/sdb/proj/Thibault/marasmius/assemblies/1_masurca/marro/CA/marro.filtered -1 /mnt/sdb/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M21-Marro_S114-S115_L003_R1_001.fastq.gz -2 /mnt/sdb/proj/Thibault/marasmius/data/illumina/files/WE-3693/231124_A00181_0697_AHVCYWDSX7/WE-3693-M21-Marro_S114-S115_L003_R2_001.fastq.gz -S marro1.sam -p 16 > marro1-bowtie2.log 2>&1 &
		#samtools
		conda activate samtools
		for i in $(ls *.sam); do out=$(echo $i | sed 's/\.sam/_sorted\.bam/'); samtools sort -@ 4 $i > $out;done
		conda activate purge_haplotigs
		for i in $(ls *.bam); do strain=$(echo $i | cut -d'_' -f1); purge_haplotigs hist -b $i  -g ../../1_hifiasm/$strain\.bp.p_ctg.gfa.filtered.fna  -t 16;done >PH_hist.log 2>&1 & #mostly only one peak, I take haplotype 1 and resume with purge_dups
	#duplications, last approach, removing small scaffolds
		#test on bamy1
		cd assemblies/5_funnaotate/busco_results/BUSCO_baemy1.hifiasm.blob.purged.ntlink.sort.masked.fa/run_agaricales
		grep Duplicated full_table.tsv | awk '{gsub("baemy1_s_","",$3); print $1, $3}' | sort -k1,1n -k2,2n | uniq | awk 'NR % 2 == 1 {first[$2]++; next}NR % 2 == 0 {last[$2]++}END {print "first";for (f in first) printf "%s %d\n", f, first[f] | "sort -k2,2nr";close("sort -k2,2nr");print "last";for (l in last) printf "%s %d\n", l, last[l] | "sort -k2,2nr";close("sort -k2,2nr")}' #for duplicated genes, select scaffold number, then for each duplication classifying the low and the high scaffold number (as proxy for size), and ranking the short and long scaffolds with the number of duplicated genes on
		cd ../../../../6_cleaned/ #>now in ../duplications
		grep ">baemy1_s_101" baemy1.fa -n #getting line number of scaff 101
		head -n 1055818 baemy1.fa > baemy1_100s.fa #cutting after scaff 100
		mkdir busco_results
		busco -i baemy1_100s.fa -m geno --offline --download_path /mnt/sdb/proj/Thibault/busco_downloads/lineages/ -l agaricales_odb10 -c 40 --out_path busco_results -q > busco_results/baemy1_100s.log 2>&1 & #C:97.6%[S:89.6%,D:8.0%],F:0.3%,M:2.1%,n:3870,E:4.3%
		cat ../../../5_funannotate/busco_results/BUSCO_baemy1.hifiasm.blob.purged.ntlink.sort.masked.fa/short_summary.specific.agaricales_odb10.BUSCO_baemy1.hifiasm.blob.purged.ntlink.sort.masked.fa.txt 		#C:97.6%[S:88.5%,D:9.1%],F:0.3%,M:2.1%,n:3870,E:4.3% >damn! gonna try with even less!
		#																																																		 C:97.0%[S:90.7%,D:6.3%],F:0.3%,M:2.7%,n:3870,E:4.3% >with 80 scaffolds, missing increase
		#																																																		 C:97.5%[S:90.2%,D:7.3%],F:0.3%,M:2.2%,n:3870,E:4.3% >I will stick to that one, looking at other stains
		#works more or less with mycal2
		#I take haplotype1 instead
#scaffolding
	#yaHS
		conda create -n yahs yahs
		conda activate samtools
		samtools faidx pr_017_003.bp.p_ctg.gfa.filtered.fna
		conda deactivate
		yahs pr_017_003.bp.p_ctg.gfa.filtered.fna ../../mapping/pr_017_003_sorted.bam #core dumped, mapping is probably deprecated?
		#mapping
			conda create -n pbmm2 pbmm2
			conda activate pbmm2
			cp pr_017_003.bp.p_ctg.gfa.filtered.fna 003.fa
		pbmm2 align 003.fa ../../data/ccsreads/pr_017_003/pr_017_003.ccsreads.fastq.gz --sort -j 12 -J 4 > 003_pbmm2.bam 2>pbmm2.log &
#remove duplicates with merqury
		cd assemblies
		mkdir merqury
		cd merqury
		sh /mnt/sdb/proj/Thibault/miniconda3/envs/merqury/share/merqury/best_k.sh 60M #get th k-mer size with 60M
		mkdir mycal1
		cd mycal1
		nohup meryl k=8 count ../../../data/ccsreads/mycal1_002/pr_017_002.ccsreads.fastq.gz output mycal1_reads.meryl > mycal1_reads_meryl.logerr 2>&1 & #get the meryl for reads, same for others
		cp ../../hifiasm/pr_017_002.bp.p_ctg.gfa.filtered.fna mycal1.fasta
		cp ../../hifiasm/pr_017_002.bp.hap1.p_ctg.gfa mycal1.hap1.gfa
		cp ../../hifiasm/pr_017_002.bp.hap2.p_ctg.gfa mycal1.hap2.gfa
		nohup meryl k=8 count mycal1.hap1.gfa output mycal1_hap1.meryl > mycal1_hap1_meryl.logerr 2>&1 &
		nohup meryl k=8 count mycal1.hap2.gfa output mycal1_hap2.meryl > mycal1_hap2_meryl.logerr 2>&1 &
		nohup merqury.sh mycal1_reads.meryl mycal1_hap1.meryl mycal1_hap2.meryl mycal1.fasta run1 > run1.logerr 2>&1 &
#install augustus (outdated)
	#https://bioinf.uni-greifswald.de/augustus/downloads/
	#option1 (done):
	#git clone https://github.com/Gaius-Augustus/Augustus
	#option 2:conda install -c bioconda augustus >can be problematic with python version
	#augustus: error while loading shared libraries: libboost_iostreams.so.1.85.0: cannot open shared object file: No such file or directory 
	#version libboost_iostreams.so.1.74.0 in miniconda3/envs/funannotate/lib
	#copying version 1.85 from busco conda
	#ln -s /mnt/sdb/proj/Thibault/miniconda3/envs/busco/lib/libboost_iostreams.so.1.85.0 /mnt/sdb/proj/Thibault/miniconda3/envs/funannotate/lib
	#augustus: error while loading shared libraries: libboost_serialization.so.1.85.0: cannot open shared object file: No such file or directory
	#same issue, same solution
	#ln -s /mnt/sdb/proj/Thibault/miniconda3/envs/busco/lib/libboost_serialization.so.1.85.0 /mnt/sdb/proj/Thibault/miniconda3/envs/funannotate/lib	
#install braker
	http://topaz.gatech.edu/GeneMark/braker.html
	#option 1 manually (very tedious)
	#option 2 (done) with conda (version 1.9)
	conda activate funannotate
	conda install braker
	#option 3 with singularity: could be installed with conda but not sure its gonna work

