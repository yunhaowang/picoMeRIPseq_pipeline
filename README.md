# In-house Python and R scripts for analyzing picoMeRIP-seq data (version 0.1)
Time-stamp: <2026-4-10 Yunhao Wang, Email: yunhaowang@126.com>


## Introduction

MeRIP-seq (N6-methyladenosine/m6A specific methylated RNA immunoprecipitation followed by high-throughput sequencing) is an experimental approach to map RNAs with specific chemical/epigenetic modifications, such as N(6)-methyladenosine (m6A). MeRipBox is a bioinformatics tool/pipeline to analyze MeRIP-seq data, and it includes 6 sections: (1) Quality control of sequencing reads; (2) Align sequencing reads to reference genome; (3) Remove multiply-aligned reads, PCR duplicates and ribosomal RNA-derived reads, and make bigWig file for data visualization; (4) Quantify gene and transcript expression; (5) Call peaks; (6) Search consensus motifs in identified peaks.


## Prerequisite

- Linux system

- Python 3.9

- Perl 5.22

- R 3.5

- SRA Toolkit (version: 2.10.9) (https://github.com/ncbi/sra-tools).

- FastQC (version: 0.12.1) (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

- Cutadapt (version: 4.8) (https://cutadapt.readthedocs.io/en/latest/guide.html).

- HISAT2 (version: 2.1.0) (http://daehwankimlab.github.io/hisat2/).

- SAMtools (version: 1.18) (http://www.htslib.org/).

- StringTie (version: 3.0.3) (https://ccb.jhu.edu/software/stringtie/).

- BEDtools (version: 2.30.0) (https://bedtools.readthedocs.io/en/latest/).

- deepTools (version: 3.5.6) (https://deeptools.readthedocs.io/en/latest/index.html).

- MACS3 (version: 3.0.3) (https://github.com/macs3-project/MACS).

- HOMER2 (version: 5.1) (http://homer.ucsd.edu/homer/). 

- MetaPlotR 46 (https://github.com/olarerin/metaPlotR). The installation path (it is “~/metaPlotR-master/” in this protocol) will be used in the metagene analysis.

- IGV (version: 2.10.3) (https://software.broadinstitute.org/software/igv/home). Install IGV on a desktop or laptop computer.

- In-house Python and R scripts (https://github.com/yunhaowang/picoMeRIPseq_pipeline/Utilities/). These scripts should be downloaded and placed in the working directory.


## Install and Run

1. Download the bioinformatics tools using the link provided above.

2. Install and define these tools in the PATH, meaning that all command lines provided in this protocol can be run directly under the working directory (commands are prefixed with a ‘$’ character).

## Pipeline

### Data preparation

1. Obtain the picoMeRIP-seq sequencing data (paired-end) using the SRA Toolkit from the NCBI BioProject under accession number PRJNA766854. The SRA IDs and corresponding sample IDs are listed in the following table. For each of the six SRA IDs, use ‘fastq-dump’ to download and convert the data into FASTQ format. The ‘fastq-dump’ command with ‘--split-3’ generates ‘_1.fastq.gz’ (read 1) and ‘_2.fastq.gz’ (read 2) for paired-end data. Below is an example using “SRR16096235”.

$ fastq-dump --split-3 --gzip SRR16096235

Table

SRA_ID Sample_ID

SRR16096235 Liver_Input

SRR16096240 Liver_IP_rep1

SRR16096241 Liver_IP_rep2

SRR16096242 Oocyte_Input

SRR24179413 Oocyte_IP_rep1

SRR24179414 Oocyte_IP_rep2

2. Rename the downloaded SRR files according the sample IDs as shown in the above table. Below is an example for the sample “Liver_Input”. 

$ mkdir ./Fastq/

$ mv SRR16096235_1.fastq.gz ./Fastq/Liver_Input_1.fastq.gz

$ mv SRR16096235_2.fastq.gz ./Fastq/Liver_Input_2.fastq.gz

3. Download the reference genome sequence file (version: mm10/GRCm38) from the UCSC Genome Browser and uncompress it.

$ wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz

$ tar -zxvf chromFa.tar.gz

$ rm chromFa.tar.gz

$ mkdir -p ./RefGenome/chromFa/

$ mv chr*.fa ./RefGenome/chromFa/

4. Concatenate all chromosome FASTQ files into a single reference genome file and generate the index (FAI format) file.

$ ls ./RefGenome/chromFa/chr*.fa | xargs cat >./RefGenome/genome_mm10.fa

$ samtools faidx ./RefGenome/genome_mm10.fa

5. Download the gene annotation library file from GENCODE. Use the version M24, which corresponds to the reference genome version mm10/GRCm38.

$ wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gtf.gz

$ gzip -d gencode.vM24.annotation.gtf.gz

$ mkdir ./Anno/

$ mv gencode.vM24.annotation.gtf ./Anno/

CRITICAL STEP: The version of the gene annotation library must match the version of reference genome.

6. Generate a BED-format gene annotation file.

$ awk '{if($3=="transcript")print $1"\t"$4"\t"$5"\t.\t.\t"$7}' ./Anno/gencode.vM24.annotation.gtf | sort -k1,1 -k2,2n >./Anno/gencode.vM24.annotation.bed

7. Extract ribosomal RNA annotations for downstream read filtering.

$ grep rRNA ./Anno/gencode.vM24.annotation.gtf | awk 'BEGIN{FS="\t"};{if($3=="gene")print $1"\t"$4"\t"$5"\t"$9"\t.\t"$7}' | sort -k1,1 -k2,2n -k3,3n >./Anno/gencode.vM24.annotation.rRNA.sort.bed


### Quality assessment of raw sequencing data

1.  Assess the quality of raw sequencing reads (FASTQ format) using FastQC. Below is an example for the sample “Liver_Input”.

$ fastqc ./Fastq/Liver_Input_1.fastq.gz

$ fastqc ./Fastq/Liver_Input_2.fastq.gz

CRITICAL STEP: Examine the FastQC reports, paying special attention to the “Overrepresented sequences” and “Adapter Content” parts. This will determine whether adapter trimming is required

2. (Optional) Trim sequencing adapters from read ends using Cutadapt. In the example dataset used in this protocol, Illumina TruSeq adapters werre detected at the 3’ end of the reads. The following command removes them for paired-end reads. Below is an example for the sample “Liver_Input”.

$ cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -j 20 -m 20 --max-n 0.01 --trim-n -o ./Fastq/Liver_Input_1.trim.fastq.gz -p ./Fastq/Liver_Input_2.trim.fastq.gz ./Fastq/Liver_Input_1.fastq.gz ./Fastq/Liver_Input_2.fastq.gz

CRITICAL STEP: The adapter sequences to be trimmed should be determined empirically for your own sequencing data, rather than copied directly from this protocol.

3. Re-run FastQC on trimmed reads to confirm quality improvement. Below is an example for the sample “Liver_Input”.

$ fastqc ./Fastq/Liver_Input_1.trim.fastq.gz

$ fastqc ./Fastq/Liver_Input_2.trim.fastq.gz


### Read alignment 

1. Build the HISAT2 index for alignment.

$ mkdir -p ./Hisat2Index/

$ hisat2-build -p 20 ./RefGenome/genome_mm10.fa ./Hisat2Index/genome_mm10

2. Align trimmed reads to the reference genome using HISAT2. Below is an example for the sample “Liver_Input”.

$ mkdir -p ./Hisat2/

$ hisat2 -p 20 --no-mixed --no-discordant -x ./Hisat2Index/genome_mm10 -1 ./Fastq/Liver_Input_1.trim.fastq.gz -2 ./Fastq/Liver_Input_2.trim.fastq.gz -S ./Hisat2/Liver_Input.sam >./Hisat2/Liver_Input.stats

CRITICAL STEP: Examine the alignment statistics (the files with the suffix “.stats”) before proceeding. If the overall alignment rate is very low (e.g., <50%), it indicates substantial contamination in the biological samples. 


### Read filtering and statistics

1. Extract uniquely-aligned reads using SAMtools. Below is an example for the sample “Liver_Input”.

$ samtools view -b -q 30 -@ 20 -o ./Hisat2/Liver_Input.uniq.bam ./Hisat2/Liver_Input.sam

$ samtools stats -@ 20 ./Hisat2/Liver_Input.uniq.bam >./Hisat2/Liver_Input.uniq.stats 

CRITICAL STEP: The files with the suffix “.uniq.stats” contain the number of read pairs that are uniquely aligned to reference genome. Note that unique alignment is determined by a mapping quality score (MAPQ ≥30), which may differ from the unique alignment rates reported by HISAT2.

2. Remove PCR duplicates using SAMtools. Below is an example for the sample “Liver_Input”.

$ samtools fixmate -@ 20 -m ./Hisat2/Liver_Input.uniq.bam ./Hisat2/Liver_Input.uniq.fm.bam

$ samtools sort -O BAM -@ 20 -o  ./Hisat2/Liver_Input.uniq.fm.sort.bam ./Hisat2/Liver_Input.uniq.fm.bam

$ samtools index -@ 20 ./Hisat2/Liver_Input.uniq.fm.sort.bam

$ samtools markdup -r -O BAM -@ 20 -s ./Hisat2/Liver_Input.uniq.fm.sort.bam ./Hisat2/Liver_Input.uniq.fm.sort.rd.bam

$ samtools index -@ 20 ./Hisat2/Liver_Input.uniq.fm.sort.rd.bam

$ samtools stats -@ 20 ./Hisat2/Liver_Input.uniq.fm.sort.rd.bam >./Hisat2/Liver_Input.uniq.fm.sort.rd.stats

CRITICAL STEP: The files with the suffix “.uniq.fm.sort.rd.stats” contain the number of read pairs after PCR duplication removal. In picoMeRIP-seq experiments, lower starting RNA amounts typically require more PCR cycles during library preparation, resulting in higher PCR duplication rates. The PCR duplication rate can therefore be used to optimize PCR cycle numbers during library preparation.

3. Remove the reads that are aligned to ribosomal RNAs using BEDtools. Below is an example for the sample “Liver_Input”.

$ bedtools intersect -a ./Hisat2/Liver_Input.uniq.fm.sort.rd.bam -b ./Anno/gencode.vM24.annotation.rRNA.sort.bed -v >./Hisat2/Liver_Input.uniq.fm.sort.rd.rr.bam

$ samtools index -@ 20 ./Hisat2/Liver_Input.uniq.fm.sort.rd.rr.bam

$ samtools stats -@ 20 ./Hisat2/Liver_Input.uniq.fm.sort.rd.rr.bam >./Hisat2/Liver_Input.uniq.fm.sort.rd.rr.stats

CRITICAL STEP: The files with the suffix “.uniq.fm.sort.rd.rr.stats” contain the number of uniquely-aligned, PCR-duplicate-free and rRNA-filtered reads. This number indicates the efficiency of the rRNA depletion step during library preparation. Importantly, these reads are the effective reads that will be used for all subsequent analyses.


### Gene transcript quantification for Input samples

1. Quantify gene and isoforme expression levels in Input samples using StringTie. Below is an example for the sample “Liver_Input”.

$ mkdir -p ./StringTie/

$ stringtie ./Hisat2/Liver_Input.uniq.fm.sort.rd.rr.bam -G ./Anno/gencode.vM24.annotation.gtf -o ./StringTie/Liver_Input.gtf -A ./StringTie/Liver_Input.gene.tab -p 20 -e

$ python py_extract_stringtie_expr.py ./StringTie/Liver_Input.gtf ./StringTie/Liver_Input.gene.tab ./StringTie/Liver_Input.gene.txt ./StringTie/Liver_Input.iso.txt


### Peak calling

1. Obtain the effective genome size from the UCSC Table Browser for use with MACS3 as follows. Open the UCSC Table Browser using the link (http://genome.ucsc.edu/cgi-bin/hgTables). Select the following parameters: “Genome” is “Mouse - Dec.2011 (GRCm38/mm10)”; “Group” is “Genes and Gene Predictions”; “Track” is “All GENCODE VM24”; “Table” is “Comprehensive (wgEncodeGencodeCompVM24)”; “Region” is “Genome”. Click the button “Summary/statistics”. A new window showing summary statistics will appear. Use the number of “block total” (231481587 in this protocol) as the effective genome size.

CRITICAL STEP: For MACS3 peak calling, the effective genome size (the parameter “-g”) should be defined based on the transcriptome size rather than the genome size.

2. Call m6A peaks using MACS3, comparing the IP sample with its corresponding Input sample. Below is an example for the sample “Liver_IP_rep1”.

$ mkdir -p ./Macs3/

$ macs3 callpeak -t ./Hisat2/Liver_IP_rep1.uniq.fm.sort.rd.rr.bam -c ./Hisat2/Liver_Input.uniq.fm.sort.rd.rr.bam -f BAM -g 231481587 --keep-dup all --outdir Macs3 -n Liver_IP_rep1 --nomodel --call-summits >./Macs3/Liver_IP_rep1.stats

CRITICAL STEP: The files with the suffix “.stats” are summaries of peak calling. For peak calling, the default cutoff is q value <0.05. The files with the suffix “_peaks.xls” are significantly-enriched peaks, which can be opened in Microsoft Excel. The files with the suffix “_peaks.narrowPeak” are genomic intervals of the significantly-enriched peaks, which can be loaded into IGV for visualization. The files with the suffix “_summits.bed” are peak summit positions of the significantly-enriched peaks, which will be used for the subsequent peak annotation, motif search, and metagene analysis.


### Peak annotation

1. Identify overlaps between peak submits and gene strands using BEDtools. Assign strand information to peak submits using an in-house Python script. Below is an example for the sample “Liver_IP_rep1”.

$ mkdir -p ./PeakAnno/

$ bedtools intersect -a ./Macs3/Liver_IP_rep1_summits.bed -b ./Anno/gencode.vM24.annotation.bed -wao | cut -f1-5,11 | uniq -c >./PeakAnno/Liver_IP_rep1.txt

$ python py_assign_strand2peak.py ./PeakAnno/Liver_IP_rep1.txt ./PeakAnno/Liver_IP_rep1.bed 

CRITICAL STEP: This is a key step because the peak submits reported by MACS3 lack strand information, which is required for peak annotation, motif search and metagene analysis. 

2. Annotate peaks using Homer2. Below is an example for the sample “Liver_IP_rep1”.

$ annotatePeaks.pl ./PeakAnno/Liver_IP_rep1.bed mm10 -cpu 20 -annStats ./PeakAnno/Liver_IP_rep1.annStats >./PeakAnno/Liver_IP_rep1.ann 


### Motif search

1. Create subfolders to store motif analysis results. Analyze motifs surrounding peak summits using HOMER2. Below is an example for the sample “Liver_IP_rep1”.

$ mkdir -p ./PeakMotif/Liver_IP_rep1/preparsed/

$ findMotifsGenome.pl ./PeakAnno/Liver_IP_rep1.bed mm10 ./PeakMotif/Liver_IP_rep1/ -len 5,6,7 -rna -p 20 -preparse -preparsedDir ./PeakMotif/Liver_IP_rep1/preparsed/
$ python py_make_RNA_motif_seqLogo.py PeakMotif/Liver_IP_rep1/homerResults/motif1.motif PeakMotif/Liver_IP_rep1.motif1.mtx PeakMotif/Liver_IP_rep1.motif1.pdf

CRITICAL STEP: If necessary, the parameter “-size” (default: 200 bp) can be adjusted based on the length of peaks called by MACS3. In most cases, 200 bp works well. Optionally, for better visulization of motifs output by HOMER2, the in-house Python code (‘py_make_RNA_motif_seqLogo.py’) can be used to re-generate PDF-format sequence logos for motifs of interest (e.g., motif1 in this protocol). 


### Metagene analysis

In this protocol, we assume the installation path of MetaPlotR tool is “~/metaPlotR-master/”.

1. Download the gene annotation file in genePred format using UCSC Table Browser. Open the UCSC Table Browser using the link (http://genome.ucsc.edu/cgi-bin/hgTables). Select the following parameters:  “Genome” is “Mouse - Dec.2011 (GRCm38/mm10)”; “Group” is “Gens and Gene Predictions”; “Track” is “All GENCODE VM24”; “Table” is “Comprehensive (wgEncodeGencodeCompVM24)”; “Region” is “Genome”; “Output format” is “All fields from selected table”. Under “Output filename”, fill “gencode.vM24.annotation.genePred”. Click the button “Get output” to download the file “gencode.vM24.annotation.genePred”. Transfer the file “gencode.vM24.annotation.genePred” to the working folder.

CRITICAL STEP: To run MetaPlotR, a gene annotation file (genePred format) is required. 

2. Prepare annotation files for running MetaPlotR.

$ mkdir ./MetaGene/

$ mv gencode.vM24.annotation.genePred ./MetaGene/

$ perl ~/metaPlotR-master/make_annot_bed.pl --genomeDir ./RefGenome/chromFa/ --genePred ./MetaGene/gencode.vM24.annotation.genePred >./MetaGene/gencode.vM24.annotation.genePred.bed

$ sort -k1,1 -k2,2n ./MetaGene/gencode.vM24.annotation.genePred.bed > ./MetaGene/gencode.vM24.annotation.genePred.sort.bed

$ perl ~/metaPlotR-master/size_of_cds_utrs.pl --annot ./MetaGene/gencode.vM24.annotation.genePred.sort.bed > ./MetaGene/region_sizes.gencode.vM24.annotation.txt

3. Identify overlaps between peak summits and gene annotations using BEDtools. Below is an example for the sample “Liver_IP_rep1”.

$ bedtools intersect -a ./PeakAnno/Liver_IP_rep1.bed -b ./MetaGene/gencode.vM24.annotation.genePred.sort.bed -sorted -wo -s >./MetaGene/Liver_IP_rep1.bed

4. Calculate relative distances using MetaPlotR. Below is an example for the sample “Liver_IP_rep1”.

$ perl ~/metaPlotR-master/rel_and_abs_dist_calc.pl --bed ./MetaGene/Liver_IP_rep1.bed --regions ./MetaGene/region_sizes.gencode.vM24.annotation.txt >./MetaGene/Liver_IP_rep1.dist.txt

5. Select the gene isoform with the highest expression level per gene using an in-house Python script. Below is an example for the sample “Liver_IP_rep1”.

$ python py_choose_highestExpr_isoform.py ./StringTie/Liver_Input.iso.txt ./MetaGene/Liver_IP_rep1.dist.txt ./MetaGene/Liver_IP_rep1.dist.highIso.txt

Only one gene isoform per gene is used for plotting metagene profile. In this protocol, we choose the isoform with the highest expression level.

6. Scale the length of the 5’ UTR and 3’ UTR relative to the CDS using an in-house Python script. Below is an example for the sample “Liver_IP_rep1”.

$ python py_scale_UTR.py 0.85 2.8 ./MetaGene/Liver_IP_rep1.dist.highIso.txt ./MetaGene/Liver_IP_rep1.dist.highIso.sc.txt

CRITICAL STEP: Before running this step, the relative length of the 5’ UTR and 3’ UTR compared to the CDS should be separately calculated based on all expressed gene isoforms. The relative length of the CDS is defined as 1.0, and its coordinates on the X-axis range from 1.0 to 2.0. The relative length of the 5’ UTR is RatioA = (average 5’ UTR length) / (average CDS length). The relative length of the 3’ UTR is RatioB = (average 3’ UTR length) / (average CDS length). The first parameter in the command (‘python py_scale_UTR.py’) is “1 – RatioA”;  the second parameter is “2 + RatioB”. In this protocol, the average lengths of CDS, 5’ UTR and 3’ UTR are 1000 bp, 150 bp and 800 bp, respectively. Thus, the first and second parameters are 0.85 and 2.8, respectively.

7. Generate a matrix for plotting metagene profiles using an in-house Python script.

$ ls ./MetaGene/*dist.highIso.sc.txt | sort >./MetaGene/metagene.fofn

$ python py_make_mtx_metagene.py ./MetaGene/metagene.fofn ./MetaGene/metagene.txt 

8. Plot the metagene profile using an in-house R script.

$ Rscript --vanilla run_plot_metagene_density.R ./MetaGene/metagene.txt ./MetaGene/metagene.pdf
 

### Genome browser visualization

The picoMeRIP-seq read density distribution (e.g., bigWig format), the peaks identified by MACS3 (e.g., BED format), and the gene annotation file (e.g., GTF format) can be loaded to the UCSC Genome Browser or IGV for visualization. 

1. Convert BAM format to bigWig format using deepTools. Below is an example for the sample “Liver_Input”.

$ bamCoverage -b ./Hisat2/Liver_Input.uniq.fm.sort.rd.rr.bam -bs 10 -p 20 --normalizeUsing RPKM -o ./Hisat2/Liver_Input.uniq.fm.sort.rd.rr.RPKM.bigWig 


## Citation

The manuscript is currently under review. The citation information will be updated once the manuscript is published. Please cite  Yunhao Wang's GitHub website at this moment.
