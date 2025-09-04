#!/usr/bin/env bash
eval "$(conda shell.bash hook)"
cpu=70
echo "cpu = $cpu"
#for i in $(ls ../../4_ntlink/hap1/ | grep -v "busco"); do
for i in $(ls ../../assemblies/6_cleaned/*fa|sed 's#.*/\([^/.]*\)\..*#\1#'); do 
	echo $i
	mkdir $i/
	echo copying input
	#cp ../../4_ntlink/hap1/$i/*fa.k32.w100.z1000.ntLink.scaffolds.gap_fill.fa $i/$i.fa
	echo activating conda funannotate
	conda activate funannotate
	#echo funannotate sort
	#funannotate sort -i $i/$i.fa -o $i/$i.sort.fa -b $i\_s --minlen 0 > $i/$i.sort.log 2>&1
	#echo funannotate mask
	#funannotate mask -i $i/$i.sort.fa -m repeatmasker -m repeatmodeler -s agaricales -l /mnt/sdb/proj/Thibault/marasmius/db/repbase_all_2022-11.fasta --cpus $cpu -o $i/$i.sort.masked.fa > $i/$i.mask.log 2>&1
	echo funannotate predict
	#funannotate predict -i $i.fa -o $i/ -s "$i" --busco_db agaricomycota --ploidy 2 --cpus $cpu \> $i/$i.predict.log 2\>\&1
	funannotate predict -i ../../assemblies/6_cleaned/$i.fa -o $i/ -s "$i" --busco_db basidiomycota --ploidy 2 --cpus $cpu > $i/$i.predict.log 2>&1
	echo funannotate iprscan
	funannotate iprscan -i $i/ -m local -c $cpu --iprscan_path /mnt/sdb/proj/Thibault/interproscan-5.73-104.0/ >$i/$i.iprscan.log 2>&1
	echo funannotate annotate
	funannotate annotate -i $i/ -o $i/ --busco_db basidiomycota --cpus $cpu >$i/$i.annotate.log 2>&1
	echo antismash
	conda activate antismash
	antismash $i/predict_results/$i.gbk -t fungi -c $cpu --genefinding-tool glimmerhmm --genefinding-gff3 $i/predict_results/$i.gff3 --output-dir $i/antismash --output-basename $i >$i/$i.antismash.log 2>&1
	conda deactivate
	mkdir $i/phobius
	echo phobius
	conda activate seqkit
	seqkit split --by-id --by-id-prefix '' $i/annotate_misc/genome.proteins.fasta -O $i/phobius/
	conda deactivate
	conda activate perl-findbin
	#for j in $(ls $i/phobius/*fasta); do perl /mnt/sdb/proj/Thibault/phobius/phobius.pl -short $j > $(echo $j | sed 's/fasta/phobius/'); done #only one protein at the teome, toos low
	ls $i/phobius/*fasta | xargs -n 1 -P $cpu -I{} sh -c 'perl /mnt/sdb/proj/Thibault/phobius/phobius.pl -short "$1" > "${1%.fasta}.phobius"' _ {} #dos it by batches of $cpu
	conda deactivate
	tail -q -n1 $i/phobius/*phobius > $i/phobius/$i.phobius
	echo annotate to combine
	funannotate annotate -i $i/ -o $i/ --antismash $i/antismash/$i.gbk --phobius $i/phobius/$i.phobius --cpus $cpu >$i/$i.combine.log 2>&1
done
