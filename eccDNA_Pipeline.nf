//scripts
test_script="$workflow.projectDir/eccDNA_Workflow/test/test.py"

proj_name=params.proj_name
study_name=params.study_name
user_name=params.user_name
output_dir=params.output_dir

sample_sheet = params.sample_sheet

trimmomatic_jar = params.trimmomatic_jar
picard_jar = params.picard_jar

refGenome_Fa = params.refGenome_Fa
refAnno_GTF  = params.refAnno_GTF
refGenome_Size = params.refGenome_Size

samples = Channel.fromPath(sample_sheet)
            .splitCsv(header:true, sep:",")
            .map {row->tuple(row.sample_id, row.sample_name,file(row.fastq1), file(row.fastq2), row.output_dir)}.view()

samples_1 = Channel.fromPath(sample_sheet)
            .splitCsv(header:true, sep:",")
            .map {row->tuple(row.sample_id, row.sample_name,file(row.fastq1), file(row.fastq2), row.output_dir)}.view()

process fastqc{
    cpus 2
    memory "10G"

    publishDir "${output_dir}/${id}/qc", mode:'copy'
    
    input:
        set val(id), val(name),file(fq1), file(fq2), val(output_dir) from samples
    output:
        set val(id), val(name),file("${fq1.baseName}_fastqc.html"),file("${fq2.baseName}_fastqc.html") into qc_samples
        file "*_fastqc.zip" into fastqcOutputZip
        file "*_fastqc.html" into fastqcOutputHtml
        file 'version.txt' into version_ch1
        file 'manifest.txt' into manifest_ch1
    
    script:
        """
            mkdir -p ./qc
            fastqc -f fastq ${fq1} -o ./
            fastqc -f fastq ${fq2} -o ./
            fastqc --version > version.txt
            ls *_fastqc.zip > manifest.txt
        """
}

process trimfastq{
    cpus 4
    memory "10G"

    input:
        set val(id), val(name),file(fq1), file(fq2), val(output_dir) from samples_1
    output:
        set val(id), val(name),file("${id}_1_tr.fq.gz"),file("${id}_2_tr.fq.gz"),val(output_dir) into trimedfastq
        file "trimmomatic.report.txt" into trimmomaticQC

    script:
        """
        module load trimmomatics/0.39

        java -jar ${trimmomatic_jar} PE -threads 4 -phred33 \
            ${fq1} \
            ${fq2} \
            ${id}_1_tr.fq.gz \
            ${id}_1_tr_unparied_fq.gz \
            ${id}_2_tr.fq.gz \
            ${id}_2_tr_unparied_fq.gz \
            ILLUMINACLIP:${params.adapter_file}:2:30:10:2:keepBothReads \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36 2> trimmomatic.report.txt
        """
}

process bwa_algnment{
    cpus 6
    memory "10G"

    //publishDir "${output_dir}/${id}/bam", mode:'copy'

    input:
        set val(id), val(name),file(fq1),file(fq2),val(output_dir) from trimedfastq
    output:
        set val(id),val(name),file("${id}_hg38_sort.bam"), file("${id}_hg38_sort.bam.bai"),val(output_dir) into sortBams
    
    script:
        """
            module load bwa samtools
            bwa mem -t 6 $refGenome_Fa -Y ${fq1} ${fq2} | samtools view -S -bS -> ${id}_hg38.bam
            samtools sort ${id}_hg38.bam -o ${id}_hg38_sort.bam -@ 6
            rm ${id}_hg38.bam
            samtools index ${id}_hg38_sort.bam

        """
}

process removedup{
    cpus 6
    memory "10G"

    publishDir "${output_dir}/${id}/bam", mode:'copy'

    input:
        set val(id), val(name),file(bam),file(bai),val(output_dir) from sortBams
    output:
        set val(id),val(name),file("${id}_hg38_sort_du.bam"), file("${id}_du_metrics.txt"), file("${id}.mk.log") into du_Bams, du_Bams_1, du_Bams_2
        file "*_du_metrics.txt" into picard_DuMatrix
    
    script:
        """
            module load picard
            java -jar ${picard_jar} MarkDuplicates I=${bam} O=${id}_hg38_sort_du.bam M=${id}_du_metrics.txt REMOVE_DUPLICATES=True 2> ${id}.mk.log
        """
}

process sotredByName{
    cpus 6
    memory "10G"

    publishDir "${output_dir}/${id}/bam", mode:'copy'

    input:
        set val(id), val(name),file(bam),file(metrics),file(log) from du_Bams
    output:
        set val(id),val(name),file("${id}_hg38_sort_du_byname.bam") into du_sortbyNameBams
    
    script:
        """
            module load samtools
            samtools view -u -f 1 -F 12 ${bam} > map_map.bam
            samtools sort -n map_map.bam -o ${id}_hg38_sort_du_byname.bam
            rm map_map.bam
        """
        //java -jar ${picard_jar} SortSam I=${bam} O=${id}_hg38_sort_du_byname.bam SORT_ORDER=queryname
}

process macs{
    cpus 2
    memory "10G"

    publishDir "${output_dir}/${id}", mode:'copy'
    input:
        set val(id),val(name),file(bam), file(metrics), file(log) from du_Bams_1
    output:
        set val(id),val(name),file("./macs2/${id}.bigwig"), file("./macs2/p_peaks.xls"), file("./macs2/p_peaks.broadPeak") into peaks
    script:
        """ 
            module load bedtools bedutil
            macs2 callpeak -t ${bam} --broad -g hs --broad-cutoff 0.1 -B -f BAMPE -np --outdir ./macs2
            cat ./macs2/p_treat_pileup.bdg |grep  -P '^chr..' > ./macs2/p_treat_pileup.bdg_clip
            sort -k1,1 -k2,2n ./macs2/p_treat_pileup.bdg_clip > ./macs2/p_treat_pileup.bdg_sort_clip
            bedGraphToBigWig ./macs2/p_treat_pileup.bdg_sort_clip ${refGenome_Size} ./macs2/${id}.bigwig
            rm -f ./macs2/p_treat_pileup.bdg_clip ./macs2/p_treat_pileup.bdg_sort_clip
        """
        //cat ./macs2/p_treat_pileup.bdg |grep  -P '^chr..' > ./macs2/p_treat_pileup.bdg_clip
        //bedtools slop -i ./macs2/p_treat_pileup.bdg -g ${refGenome_Size} -b 0 |bedClip stdin ${refGenome_Size} ./macs2/p_treat_pileup.bdg_clip
}

process getPicardMatrics{
    maxRetries 4
    errorStrategy {sleep(Math.pow(2, task.attempt)*200 as long); task.attempt < maxRetries ? 'retry':'ignore'}

    publishDir "${output_dir}/${id}/qc", mode:'copy'
    input:
        set val(id),val(name),file(bam), file(metrics), file(log) from du_Bams_2
    output:
        file "*.alignment_summary_metrics" into picardOut, picardOut_forMultiQC
    script:
        """
            module load picard
            picard CollectMultipleMetrics I=${bam} O=${bam.baseName} R=${refGenome_Fa}
        """
}

process peakanno{
    cpus 2
    memory "10G"
    publishDir "${output_dir}/${id}/macs2", mode:'copy'
    input:
        set val(id),val(name),file(bigwig), file(peaks), file(broadPeak) from peaks
    output:
        set val(id),val(name),file("p_peaks.broadPeak_anno") into anno_peaks
    script:
        """
            module load homer
            annotatePeaks.pl ${broadPeak} hg38 > p_peaks.broadPeak_anno
        """
}

process runMultiqc{
    maxRetries 2
    errorStrategy {sleep(Math.pow(2, task.attempt)*200 as long); task.attempt < maxRetries ? 'retry':'ignore'}
    
    publishDir "${output_dir}/multiqc", mode:'copy'

    input:
        file(fqc) from fastqcOutputHtml.collect()
        file(tqc) from trimmomaticQC.collect()
        file(pqc) from picardOut_forMultiQC.collect()
        file(dqc) from picard_DuMatrix.collect()
    output:
        file "multiqc_report.html" into multiQcOutReport

    script:
        """
            module load multiqc
            multiqc ./
        """
}
