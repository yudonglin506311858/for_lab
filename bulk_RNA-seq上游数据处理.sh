cd /data1/yudonglin_data/lxh/mating0330014/X101SC21115713-Z01-J011-B11-16/rawdata
mkdir clean
mkdir hisat
mkdir count
for i in {"10-炎红_FRAS220042664-2r","14-炎红_FRAS220042668-1r","7-炎红_FRAS220042661-2r","10-炎红黑_FRAS220042673-1r","14-非炎红_FRAS220042675-1r",8"-炎红_FRAS220042662-2r","11-炎红_FRAS220042665-2r","15-炎红_FRAS220042669-1r","9-炎红_FRAS220042663-2r","12-炎红_FRAS220042666-1r","15-非炎红_FRAS220042676-1r","9-炎红黑_FRAS220042672-1r","13-炎红_FRAS220042667-1r","3-炎红_FRAS220042670-1r","13-非炎红_FRAS220042674-1r","5-炎红_FRAS220042671-1r"};
do 
cutadapt --pair-filter=any --minimum-length 15 --max-n 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o ./clean/${i}_rmadp_1.fq.gz -p ./clean/${i}_rmadp_2.fq.gz ${i}_1.fq.gz ${i}_2.fq.gz >>filter.txt 2>&1
java -jar /data/yudonglin/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 40 -phred33 ./clean/${i}_rmadp_1.fq.gz ./clean/${i}_rmadp_2.fq.gz -baseout ./clean/${i}_fliter.fq.gz  AVGQUAL:20 SLIDINGWINDOW:4:15 MINLEN:15 1>>filter.txt 2>&1;
hisat2 --threads 35 -x /data/yudonglin/reference/mm10/genome -1 ./clean/${i}_fliter_1P.fq.gz -2 ./clean/${i}_fliter_2P.fq.gz -S ./hisat/${i}.sam;
samtools view -S ./hisat/${i}.sam -b > ./hisat/${i}.bam;
rm ./hisat/${i}.sam;
samtools sort ./hisat/${i}.bam -o ./hisat/${i}_sorted.bam;
samtools index ./hisat/${i}_sorted.bam;
featureCounts -T 30 -t exon -g gene_id -a /data/yudonglin/reference/Mus_musculus.GRCm38.102.gtf -o ./count/${i}.count ./hisat/${i}_sorted.bam >>~/count.txt 2>&1
conda activate pyscenic
bamCoverage --bam ./hisat/${i}_sorted.bam -o ./hisat/${i}_sorted.bam.bw  --binSize 10 -p 40
conda deactivate
#rm ./hisat/${i}.bam;
done



cd /data1/yudonglin_data/lxh/mating0330014/Result-P101SC18031367-01-J002-B2-21/1.rawdata/data
mkdir clean
mkdir hisat
mkdir count
for i in {"EN13_FKDL220022167-1a","EN15_FKDL220022169-1a","EN8_FKDL220022165-1a","EN14_FKDL220022168-1a","EN7_FKDL220022164-1a","EN9_FKDL220022166-1a"};
do 
cutadapt --pair-filter=any --minimum-length 15 --max-n 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o ./clean/${i}_rmadp_1.fq.gz -p ./clean/${i}_rmadp_2.fq.gz ${i}_1.fq.gz ${i}_2.fq.gz >>filter.txt 2>&1
java -jar /data/yudonglin/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 40 -phred33 ./clean/${i}_rmadp_1.fq.gz ./clean/${i}_rmadp_2.fq.gz -baseout ./clean/${i}_fliter.fq.gz  AVGQUAL:20 SLIDINGWINDOW:4:15 MINLEN:15 1>>filter.txt 2>&1;
hisat2 --threads 35 -x /data/yudonglin/reference/mm10/genome -1 ./clean/${i}_fliter_1P.fq.gz -2 ./clean/${i}_fliter_2P.fq.gz -S ./hisat/${i}.sam;
samtools view -S ./hisat/${i}.sam -b > ./hisat/${i}.bam;
rm ./hisat/${i}.sam;
samtools sort ./hisat/${i}.bam -o ./hisat/${i}_sorted.bam;
samtools index ./hisat/${i}_sorted.bam;
featureCounts -T 30 -t exon -g gene_id -a /data/yudonglin/reference/Mus_musculus.GRCm38.102.gtf -o ./count/${i}.count ./hisat/${i}_sorted.bam >>~/count.txt 2>&1
conda activate pyscenic
bamCoverage --bam ./hisat/${i}_sorted.bam -o ./hisat/${i}_sorted.bam.bw  --binSize 10 -p 40
conda deactivate
#rm ./hisat/${i}.bam;
done




cd /data1/yudonglin_data/yr/rawdata
mkdir clean
mkdir hisat
mkdir count
for i in {"11EB_FRAL220024682-1a","1EB_FRAL220024676-1a","4EB_FRAL220024679-1a","12EB_FRAL220024683-1a","3EB_FRAL220024678-1a","5EB_FRAL220024680-1a"};
do 
cutadapt --pair-filter=any --minimum-length 15 --max-n 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o ./clean/${i}_rmadp_1.fq.gz -p ./clean/${i}_rmadp_2.fq.gz ${i}_1.fq.gz ${i}_2.fq.gz >>filter.txt 2>&1
java -jar /data/yudonglin/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 40 -phred33 ./clean/${i}_rmadp_1.fq.gz ./clean/${i}_rmadp_2.fq.gz -baseout ./clean/${i}_fliter.fq.gz  AVGQUAL:20 SLIDINGWINDOW:4:15 MINLEN:15 1>>filter.txt 2>&1;
hisat2 --threads 35 -x /data/yudonglin/reference/mm10/genome -1 ./clean/${i}_fliter_1P.fq.gz -2 ./clean/${i}_fliter_2P.fq.gz -S ./hisat/${i}.sam;
samtools view -S ./hisat/${i}.sam -b > ./hisat/${i}.bam;
rm ./hisat/${i}.sam;
samtools sort ./hisat/${i}.bam -o ./hisat/${i}_sorted.bam;
samtools index ./hisat/${i}_sorted.bam;
featureCounts -T 30 -t exon -g gene_id -a /data/yudonglin/reference/Mus_musculus.GRCm38.102.gtf -o ./count/${i}.count ./hisat/${i}_sorted.bam >>~/count.txt 2>&1
conda activate pyscenic
bamCoverage --bam ./hisat/${i}_sorted.bam -o ./hisat/${i}_sorted.bam.bw  --binSize 10 -p 40
conda deactivate
#rm ./hisat/${i}.bam;
done
