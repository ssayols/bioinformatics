// Map as described in the PRO-seq methods section:
//
// The PRO-seq reads were adapter-clipped using cutadapt78 and trimmed and filtered
// to 15–36 bp with fastx (http://hannonlab.cshl.edu/fastx_toolkit/). Reads that
// did not map to ribosomal RNA genes were aligned to hg19 using Bowtie75 selecting
// only uniquely mapping reads with up to two mismatches. The complete raw data for
// PRO-seq in human K562 cells is available at GEO database
// (http://www.ncbi.nlm.nih.gov/geo) under accession: GSE89230.
PROJECT  = "/fsimb/groups/imb-bioinfocf/projects/sfb/roukos/sfb_roukos_2019_01_gothe_proseq/proseq"
THREADS  = "12"
RRNA_REF = "/fsimb/common/genomes/homo_sapiens/rrna/index/bowtie/1.1.2/human_all_rRNA"
HG19_REF = "/fsimb/groups/imb-bioinfocf/common-data/igenomes_reference/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome"

config {
    executor="slurm"
    queue="bcfshort"
    commands {
        cutadapt {
            walltime="02:00:00"
            procs="16"
            memory="4"
        }
        bowtie {
            walltime="02:00:00"
            procs="16"
            memory="16"
        }
        star {
            walltime="02:00:00"
            procs="16"
            memory="32"
        }
    }
}

trim = {
    doc "adapter-clip using cutadapt and trim and filter to m–M bp (avoid using fastx as the paper suggests)"

    output.dir = PROJECT + "/rawdata"
    def CUTADAPT_FLAGS =
        " --trim-n" + 
        " -b " + adapter +
        " -j " + threads +
        " -m " + min_len +
        " -M " + max_len

    transform (".fastq.gz") to (".trimmed.fastq.gz") {
        exec """
            ml pigz cutadapt/2.4 &&
            cutadapt $CUTADAPT_FLAGS -o $output $input
        """, "cutadapt"
    }
}

bowtie_rrna = {
    doc "align with bowtie (best mode), and discard everything that maps on the reference"
    output.dir = PROJECT + "/rawdata"

    def BOWTIE_FLAGS = "-q --sam --phred33-quals --chunkmbs 256" +
        " -v " + max_mismatches   +
        " -p " + threads

    def SAMTOOLS_VIEW_FLAGS = "-bhSu "
    def SAMTOOLS_SORT_FLAGS = "-O bam -@ " + threads

    transform(".fastq.gz") to (".rrna_free.fastq.gz") {
        exec """
            ml pigz bowtie/1.2.2 &&
            export TMP=/jobdir/\${SLURM_JOBID}/ &&
            TMP_FILE=\${TMP}/\$(basename $output.prefix) &&
            zcat $input | bowtie $BOWTIE_FLAGS --un \${TMP_FILE} $ref - > /dev/null &&
            pigz -cp $threads \${TMP_FILE} > $output &&
            rm \${TMP_FILE}
        """, "bowtie"
    }
}

star = {
    doc "align with star (disabling spliced alignment), selecting only uniquely mapping reads"
    output.dir = PROJECT + "/mapped"

    def STAR_FLAGS = "--runMode alignReads "        +
                     "--genomeLoad NoSharedMemory " +
                     "--outStd SAM "                +
                     "--outSAMattributes Standard " +
                     "--outSJfilterReads Unique "   +
                     "--readFilesCommand zcat "     +
                     "--outSAMunmapped None "       +
//                   "--scoreDelOpen -10000 --scoreInsOpen -10000 " +   // disable gaps and indels
                     "--alignIntronMax 1 --alignIntronMin 2 " +   // disable splicing
                     "--outFilterMismatchNmax " + max_mismatches +
                     "--runThreadN " + threads +
                     "--genomeDir " + ref +
                     "--sjdbGTFfile " + gtf

    def SAMTOOLS_VIEW_FLAGS = "-bhSu -F 256 -@ " + threads   // remove secondary alignments

    transform(".fastq.gz") to (".bam") {
        exec """
            ml star/2.6.1b_debian9 samtools/1.9 &&
            export TMP=/jobdir/\${SLURM_JOBID}/ &&
            STAR $STAR_FLAGS --outFileNamePrefix \${TMP}/STAR_output/\$(basename $input) --outTmpDir \${TMP} --readFilesIn $input |
              samtools view $SAMTOOLS_VIEW_FLAGS - |
              samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/$(basename $input)_sort - > $output &&
            mv \${TMP}/STAR_output/* $output.dir
        """, "star"
    }
}

bam2bw = {
    output.dir = PROJECT + "/tracks"

    transform(".bam") to (".scaled.bw") {
        exec """
            ml bedtools/2.27.1_debian9 samtools/1.9 kentUtils/v365 &&
            export TMP=/jobdir/\${SLURM_JOBID}/ &&

            CHRSIZES=${TMP}/\$(basename ${input.prefix}).bam2bw.chrsizes &&
            samtools idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&
            TOTAL_MAPPED=\$( samtools flagstat $input | head -n1| cut -f1 -d" ") &&
            SCALE=\$(echo "1000000/\$TOTAL_MAPPED" | bc -l) &&
            genomeCoverageBed -bg -split -scale \${SCALE} -ibam ${input} | sort -k1,1 -k2,2n > ${output.prefix}.bedgraph &&
            bedGraphToBigWig ${output.prefix}.bedgraph \${CHRSIZES} $output &&
            rm \${CHRSIZES} ${output.prefix}.bedgraph
        """
    }
}

Bpipe.run {
    "%.fastq.gz" * [
//        trim.using(adapter="TGGAATTCTCGGGTGCCAAGG", min_len=15, max_len=36, threads=THREADS) +
        trim.using(adapter="TGGAATTCTCGGGTGCCAAGG", min_len=15, max_len=100, threads=THREADS) +
        bowtie_rrna.using(max_mismatches=2, threads=THREADS, ref=RRNA_REF) +
        star.using(max_mismatches=2, threads=THREADS, ref=HG19_REF) +
        bam2bw
    ]
}

