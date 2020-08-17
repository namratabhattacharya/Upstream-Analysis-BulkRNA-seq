#!/bin/bash
echo -n "Enter genome type (human/mouse): "
read gentype
ncore=8
echo -n "Enter number samples: "
read num 
echo -n "Enter read types (PE/SE): "
read readtype
echo -e "filename\tlabel"> filelabel.txt
echo "Enter Accession number"
for (( i=1; i<$num+1; i++ ));
do
  echo "File" $i ":" 
  read sra[$i]
  echo "Label" $i "(tumor/nontumor):"
  read label
  echo -e ${sra[i]}".genes.results\t"$label >> filelabel.txt
done



if [ $gentype = "human" ]
then
  if [ ! -d "./hg" ]
  then
     mkdir hg
     echo "Import Gencode GRCh38..."
     wget -P ./hg ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.primary_assembly.genome.fa.gz
     wget -P ./hg ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.primary_assembly.annotation.gtf.gz
     gunzip -r ./hg 
  fi 
elif [ $gentype = "mouse" ]
then
  if [ ! -d "./mm" ] 
  then   
     mkdir mm
     echo "Import Gencode GRCm38..."
     wget -P ./mm ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
     wget -P ./mm ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz
     gunzip -r ./mm
  fi
else
  echo "None of the condition met"
  exit
fi

## Make STAR index
echo "Generating STAR index..."
if [ $gentype = "human" ]
then
  if [ ! -d "./sidx_hg" ]
  then
      mkdir sidx_hg
      STAR --runThreadN $ncore --runMode genomeGenerate --genomeDir sidx_hg --genomeFastaFiles ./hg/GRCh38.primary_assembly.genome.fa --sjdbGTFfile ./hg/gencode.v33.primary_assembly.annotation.gtf 
  fi
elif [ $gentype = "mouse" ]
then
  if [ ! -d "./sidx_mm" ]
  then
      mkdir sidx_mm
      STAR --runThreadN $ncore --runMode genomeGenerate --genomeDir sidx_mm --genomeFastaFiles ./mm/GRCm38.primary_assembly.genome.fa --sjdbGTFfile ./mm/gencode.vM25.primary_assembly.annotation.gtf  
  fi
else
  echo "Error: Couldn't generate STAR index..."
  exit
fi


## RSEM reference preparation
echo "Preparing RSEM reference..."
if [ $gentype = "human" ]
then
  if [ ! -d "./rsem_hg" ]
  then
     mkdir rsem_hg	
     cd rsem_hg/
     rsem-prepare-reference --gtf ../hg/gencode.v33.primary_assembly.annotation.gtf -p $ncore ../hg/GRCh38.primary_assembly.genome.fa GRCh38
     cd ..
  fi
elif [ $gentype = "mouse" ]
then
  if [ ! -d "./rsem_mm" ]
  then
     mkdir rsem_mm	
     cd rsem_mm/
     rsem-prepare-reference --gtf ../mm/gencode.vM25.primary_assembly.annotation.gtf -p $ncore ../mm/GRCm38.primary_assembly.genome.fa  GRCm38
     cd ..
  fi
else
   echo "Error: Preparaing RSEM reference..."
   exit
fi

echo "Input Samples:"
for i in "${sra[@]}"
do 
    echo $i  
done


if [ $readtype = "PE" ] 
then
  for read in "${sra[@]}"
  do
    echo "Importing" $read "..." 
    fastq-dump -I --split-files $read
    echo "Quality Checking..."
    ~/TrimGalore-0.6.5/trim_galore --fastqc --illumina --length 35 --fastqc --cores $ncore  --paired ${read}_1.fastq ${read}_2.fastq 
    if [ $gentype = "human" ]
    then
    	echo "Mapping..."
    	STAR --genomeDir ./sidx_hg --runMode alignReads --runThreadN $ncore --readFilesIn ${read}_1_val_1.fq ${read}_2_val_2.fq --quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within  --outFilterMultimapNmax 10 --outFileNamePrefix  ${read}
    	echo "Read quantification"
    	rsem-calculate-expression --paired-end --alignments --bam --no-bam-output -p 10 ${read}Aligned.toTranscriptome.out.bam  ./rsem_hg/GRCh38 ${read}
    elif [ $gentype = "mouse" ]
    then
    	STAR --genomeDir ./sidx_mm --runMode alignReads --runThreadN $ncore --readFilesIn ${read}_1_val_1.fq ${read}_2_val_2.fq --quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within  --outFilterMultimapNmax 10 --outFileNamePrefix  ${read}
    	echo "Read quantification"
    	rsem-calculate-expression --paired-end --alignments --bam --no-bam-output -p $ncore ${read}Aligned.toTranscriptome.out.bam  ./rsem_mm/GRCm38 ${read}
    else
       echo "Error: Fault in processing"  ${read} "..."
       exit
    fi
  done 
elif [ $readtype = "SE" ]
then
  for read in "${sra[@]}"
  do
    echo "Importing" $read "..." 
    fastq-dump $read
    echo "Quality Checking..."
    ~/TrimGalore-0.6.5/trim_galore --fastqc --illumina --length 35 --fastqc --cores $ncore ${read}.fastq 
    if [ $gentype = "human" ]
    then
    	echo "Mapping..."
    	STAR --genomeDir ./sidx_hg --runMode alignReads --runThreadN $ncore --readFilesIn ${read}_trimmed.fq --quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within  --outFilterMultimapNmax 10 --outFileNamePrefix  ${read}
    	echo "Read quantification"
    	rsem-calculate-expression --alignments --bam --no-bam-output -p 10 ${read}Aligned.toTranscriptome.out.bam  ./rsem_hg/GRCh38 ${read}
    elif [ $gentype = "mouse" ]
    then
    	STAR --genomeDir ./sidx_mm --runMode alignReads --runThreadN $ncore --readFilesIn ${read}_trimmed.fq --quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within  --outFilterMultimapNmax 10 --outFileNamePrefix  ${read}
    	echo "Read quantification"
    	rsem-calculate-expression  --alignments --bam --no-bam-output -p $ncore ${read}Aligned.toTranscriptome.out.bam  ./rsem_mm/GRCm38 ${read}
    else
       echo "Error: Fault in processing"  ${read} "..."
       exit
    fi
  done 
else
   echo "Error: Fault in read quantification..."
   exit
fi

mkdir Result_ReadQuant
mv *.results Result_ReadQuant/
mv filelabel.txt Result_ReadQuant/

mkdir Result_AlignPct
mv *.html  Result_AlignPct/
mv *.zip  Result_AlignPct/
if [ $readtype = "PE" ] 
then
    mv *_1_val_1.fq  Result_AlignPct/
    mv *_2_val_2.fq  Result_AlignPct/
else
    mv *_trimmed.fq  Result_AlignPct/
fi
rm -r SRR*
