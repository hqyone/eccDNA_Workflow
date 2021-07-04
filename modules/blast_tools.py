# Run BLAST for a bam file
import os
import subprocess

def getAlignmentCMD(cmdpath, db_fasta, query_fasta, eval, idcutoff, blast_out_file, hit_number=10):
    cmd = ""
    output_format_str = '-outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send qseq sseq qlen slen evalue"'
    if os.path.isfile(db_fasta) and os.path.isfile(query_fasta):
        cmd = cmdpath + " -task blastn -db " + db_fasta + " -query " + \
            query_fasta +" -perc_identity " +str(idcutoff)+ " -evalue " + str(eval) + " -num_threads 8 " + output_format_str + \
                "  -num_alignments " + \
                str(hit_number) + " > " + blast_out_file
    return cmd

def RunBLASTN(blastn, mkdb, id, db_fasta, query_fasta, blast_out_dir, eval=0.01,idcutoff=98, hit_number=5):
    try:
        print("Begin process "+id+" ...")
        blast_out_file = blast_out_dir + "/"+id+"_tRNA_blast_out.tab"
        if not os.path.isfile(db_fasta+".nhr"):
            p = subprocess.Popen("makeblastdb -in "+db_fasta+" -dbtype nucl -title "+id, shell=True, stdout=subprocess.PIPE)
            p.wait()
        blastn_cmd = getAlignmentCMD(blastn, db_fasta,query_fasta,eval,idcutoff,blast_out_file,hit_number=hit_number)
        print(blastn_cmd)
        process = subprocess.Popen(blastn_cmd, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print("Finish processing " + id)
        return blast_out_file
    except():
        return ""

