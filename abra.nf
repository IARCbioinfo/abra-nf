#!/usr/bin/env nextflow

// requires (in path):
// java
// bedtools
// bwa
// sambamba
// samtools
// abra jar

params.help = null
params.tumor_bam_folder = null
params.normal_bam_folder = null
params.bam_folder = null
params.bed = null
params.ref = null
params.abra_path = null
params.read_length = null
params.abra2 = "false"

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info '  NEXTFLOW for abra             '
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run abra_TN_pairs.nf --tumor_bam_folder tumor_BAM/ --normal_bam_folder normal_BAM/ --bed mybedfile.bed --ref ref.fasta' 
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '   When using Tumor/Normal pairs:'    
    log.info '    --tumor_bam_folder   FOLDER                  Folder containing tumor BAM files.'
    log.info '    --normal_bam_folder  FOLDER                  Folder containing matched normal BAM files.'
    log.info '   In other cases:'     
    log.info '    --bam_folder         FOLDER                  Folder containing BAM files.'
    log.info '   In all cases:'         
    log.info '    --bed                FILE                    Bed file containing intervals.'
    log.info '    --ref                FILE (with index)       Reference fasta file indexed by bwa.'
    log.info '    --abra_path          FILE                    abra.jar explicit path.'
    log.info '    --read_length        INT                     Read length (e.g.: 100).'   
    log.info 'Optional arguments:'
    log.info '   When using Tumor/Normal pairs:'  
    log.info '    --suffix_tumor       STRING                  Suffix identifying tumor bam (default: "_T").'
    log.info '    --suffix_normal      STRING                  Suffix identifying normal bam (default: "_N").'
    log.info '   In all cases:'             
    log.info '    --mem                INTEGER                 RAM used (in GB, default: 16)'
    log.info '    --threads            INTEGER                 Number of threads (default: 4)'
    log.info '    --out_folder         FOLDER                  Output folder (default: abra_BAM).'
    log.info ''
    exit 1
}

assert (params.ref != true) && (params.ref != null) : "please specify --ref option (--ref reference.fasta(.gz))"

if(params.bam_folder) {
    assert (params.bam_folder != true) && (params.bam_folder != null) : "please specify --bam_folder option (--bam_folder bamfolder)"
} else {
    assert (params.normal_bam_folder != true) && (params.normal_bam_folder != null) : "please specify --normal_bam_folder option (--normal_bam_folder bamfolder)"
    assert (params.tumor_bam_folder != true) && (params.tumor_bam_folder != null) : "please specify --tumor_bam_folder option (--tumor_bam_folder bamfolder)"
}

assert (params.bed != true) && (params.bed != null) : "please specify --bed option (--bed regions.bed)"
assert (params.abra_path != true) && (params.abra_path != null) : "please specify --abra_path option (--abra_path /path/to/abra.jar)"
assert (params.read_length != true) && (params.read_length != null) : "please specify --read_length option (--read_length 100)"


fasta_ref = file(params.ref)
fasta_ref_fai = file( params.ref+'.fai' )
fasta_ref_gzi = file( params.ref+'.gzi' )
fasta_ref_sa = file( params.ref+'.sa' )
fasta_ref_bwt = file( params.ref+'.bwt' )
fasta_ref_ann = file( params.ref+'.ann' )
fasta_ref_amb = file( params.ref+'.amb' )
fasta_ref_pac = file( params.ref+'.pac' )

bed = file(params.bed)

params.suffix_tumor = "_T"
params.suffix_normal = "_N"
params.mem = 16
params.threads = 4
params.out_folder = "abra_BAM"

try { assert file(params.bed).exists() : "\n WARNING : input bed file not located in execution directory" } catch (AssertionError e) { println e.getMessage() }

try { assert fasta_ref.exists() : "\n WARNING : fasta reference not located in execution directory. Make sure reference index is in the same folder as fasta reference" } catch (AssertionError e) { println e.getMessage() }
if (fasta_ref.exists()) {assert fasta_ref_fai.exists() : "input fasta reference does not seem to have a .fai index (use samtools faidx)"}
if (fasta_ref.exists()) {assert fasta_ref_sa.exists() : "input fasta reference does not seem to have a .sa index (use bwa index)"}
if (fasta_ref.exists()) {assert fasta_ref_bwt.exists() : "input fasta reference does not seem to have a .bwt index (use bwa index)"}
if (fasta_ref.exists()) {assert fasta_ref_ann.exists() : "input fasta reference does not seem to have a .ann index (use bwa index)"}
if (fasta_ref.exists()) {assert fasta_ref_amb.exists() : "input fasta reference does not seem to have a .amb index (use bwa index)"}
if (fasta_ref.exists()) {assert fasta_ref_pac.exists() : "input fasta reference does not seem to have a .pac index (use bwa index)"}

if (fasta_ref.exists() && params.ref.tokenize('.')[-1] == 'gz') {assert fasta_ref_gzi.exists() : "input gz fasta reference does not seem to have a .gzi index (use samtools faidx)"}

assert (params.read_length > 0) : "read length must be higher than 0 (--read_length)"

process bed_kmer_size {

    cpus params.threads

    input:
    file bed

    output:
    file "kmer_size_abra.bed" into bed_kmer

    shell:
    '''
    grep -v '^track' !{bed} | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk '{print $1"\t"$2"\t"$3}' > tmp_merged_sorted.bed
    java -Xmx4G -cp !{params.abra_path} abra.KmerSizeEvaluator !{params.read_length} /appli57/reference/ucsc.hg19.fasta kmer_size_abra.bed !{params.threads} tmp_merged_sorted.bed
    '''
}

if(params.bam_folder) {

    try { assert file(params.bam_folder).exists() : "\n WARNING : input BAM folder not located in execution directory" } catch (AssertionError e) { println e.getMessage() }
    assert file(params.bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0 : "BAM folder contains no BAM"

    // recovering of bam files
    bams = Channel.fromPath( params.bam_folder+'/*.bam' )
              .ifEmpty { error "Cannot find any bam file in: ${params.bam_folder}" }
              .map {  path -> [ path.name.replace(".bam",""), path ] }

    // recovering of bai files
    bais = Channel.fromPath( params.bam_folder+'/*.bam.bai' )
              .ifEmpty { error "Cannot find any bai file in: ${params.bam_folder}" }
              .map {  path -> [ path.name.replace(".bam.bai",""), path ] }

    // building bam-bai pairs
    bam_bai = bams
              .phase(bais)
              .map { bam, bai -> [ bam[1], bai[1] ] }
              
    process abra {

        cpus params.threads
        memory params.mem+'GB' 

        tag { bam_tag }
        
        publishDir params.out_folder, mode: 'move', pattern: '*_SV.txt'

        input:
        file bam_bai
        file bed_kmer from bed_kmer.first()
        file fasta_ref
        file fasta_ref_fai
        file fasta_ref_gzi
        file fasta_ref_sa 
        file fasta_ref_bwt
        file fasta_ref_ann
        file fasta_ref_amb
        file fasta_ref_pac

        output:
        file("${bam_tag}_abra.bam") into bam_abra
        file("${bam_tag}_SV.txt") optional true into SV_output

        shell:
        bam_tag = bam_bai[0].baseName
	if(params.abra2=="false") abraoptions="--working abra_tmp --sv tmp_SV.txt"
	else abraoptions="--tmpdir ."
        '''
        java -Xmx!{params.mem}g -jar !{params.abra_path} --in !{bam_tag}.bam --out "!{bam_tag}_abra.bam" --ref !{fasta_ref} --target-kmers !{bed_kmer} --threads !{params.threads} !{abraoptions} > !{bam_tag}_abra.log 2>&1
        if [ -f tmp_SV.txt ]; then
		mv tmp_SV.txt !{bam_tag}_SV.txt
	fi
	'''
    }


} else {

    try { assert file(params.tumor_bam_folder).exists() : "\n WARNING : input tumor BAM folder not located in execution directory" } catch (AssertionError e) { println e.getMessage() }
    assert file(params.tumor_bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0 : "tumor BAM folder contains no BAM"
    try { assert file(params.normal_bam_folder).exists() : "\n WARNING : input normal BAM folder not located in execution directory" } catch (AssertionError e) { println e.getMessage() }
    assert file(params.normal_bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0 : "normal BAM folder contains no BAM"

    // FOR TUMOR 
    // recovering of bam files
    tumor_bams = Channel.fromPath( params.tumor_bam_folder+'/*'+params.suffix_tumor+'.bam' )
              .ifEmpty { error "Cannot find any bam file in: ${params.tumor_bam_folder}" }
              .map {  path -> [ path.name.replace("${params.suffix_tumor}.bam",""), path ] }

    // recovering of bai files
    tumor_bais = Channel.fromPath( params.tumor_bam_folder+'/*'+params.suffix_tumor+'.bam.bai' )
              .ifEmpty { error "Cannot find any bai file in: ${params.tumor_bam_folder}" }
              .map {  path -> [ path.name.replace("${params.suffix_tumor}.bam.bai",""), path ] }

    // building bam-bai pairs
    tumor_bam_bai = tumor_bams
              .phase(tumor_bais)
              .map { tumor_bam, tumor_bai -> [ tumor_bam[0], tumor_bam[1], tumor_bai[1] ] }

    // FOR NORMAL 
    // recovering of bam files
    normal_bams = Channel.fromPath( params.normal_bam_folder+'/*'+params.suffix_normal+'.bam' )
              .ifEmpty { error "Cannot find any bam file in: ${params.normal_bam_folder}" }
              .map {  path -> [ path.name.replace("${params.suffix_normal}.bam",""), path ] }

    // recovering of bai files
    normal_bais = Channel.fromPath( params.normal_bam_folder+'/*'+params.suffix_normal+'.bam.bai' )
              .ifEmpty { error "Cannot find any bai file in: ${params.normal_bam_folder}" }
              .map {  path -> [ path.name.replace("${params.suffix_normal}.bam.bai",""), path ] }

    // building bam-bai pairs
    normal_bam_bai = normal_bams
              .phase(normal_bais)
              .map { normal_bam, normal_bai -> [ normal_bam[0], normal_bam[1], normal_bai[1] ] }

    // building 4-uplets corresponding to {tumor_bam, tumor_bai, normal_bam, normal_bai}
    tn_bambai = tumor_bam_bai
          .phase(normal_bam_bai)
          .map {tumor_bb, normal_bb -> [ tumor_bb[1], tumor_bb[2], normal_bb[1], normal_bb[2] ] }    
    // here each element X of tn_bambai channel is a 4-uplet. X[0] is the tumor bam, X[1] the tumor bai, X[2] the normal bam and X[3] the normal bai.


    process abra_TN {

        cpus params.threads
        memory params.mem+'GB' 

        tag { tumor_normal_tag }

        publishDir params.out_folder, mode: 'move', pattern: '*_SV.txt'

        input:
        file tn from tn_bambai
        file bed_kmer from bed_kmer.first()
        file fasta_ref
        file fasta_ref_fai
        file fasta_ref_gzi
        file fasta_ref_sa 
        file fasta_ref_bwt
        file fasta_ref_ann
        file fasta_ref_amb
        file fasta_ref_pac


        output:
    //    file("${tumor_normal_tag}${params.suffix_normal}_abra.bam") into normal_output
    //    file("${tumor_normal_tag}${params.suffix_tumor}_abra.bam") into tumor_output
        file '*_abra.bam' into bam_abra  mode flatten
        file("${tumor_normal_tag}_SV.txt") optional true into SV_output

        shell:
        tumor_normal_tag = tn[0].baseName.replace(params.suffix_tumor,"")
	if(params.abra2=="false") abraoptions="--working abra_tmp --sv tmp_SV.txt"
        else abraoptions="--tmpdir ."
               
	'''
        java -Xmx!{params.mem}g -jar !{params.abra_path} --in !{tumor_normal_tag}!{params.suffix_normal}.bam,!{tumor_normal_tag}!{params.suffix_tumor}.bam --out "!{tumor_normal_tag}!{params.suffix_normal}_abra.bam","!{tumor_normal_tag}!{params.suffix_tumor}_abra.bam" --ref !{fasta_ref} --target-kmers !{bed_kmer} --threads !{params.threads} !{abraoptions}> !{tumor_normal_tag}_abra.log 2>&1
        if [ -f tmp_SV.txt ]; then
		mv tmp_SV.txt !{tumor_normal_tag}_SV.txt
	fi
	'''
    }
}

process fixmate_sort_index {

    cpus params.threads
    memory params.mem+'GB' 

    tag { bam_tag }

    publishDir params.out_folder, mode: 'move'

    input:
    file bam_abra

    output:
    file '*abra_sorted_fixmate.bam*' into final_bam 

    shell:
    bam_tag = bam_abra.baseName
    half_mem = params.mem.intdiv(2)
    half_threads = params.threads.intdiv(2) - 1
    '''
    set -o pipefail
    sambamba sort -t !{half_threads} -m !{half_mem}G -n --tmpdir=sort_tmp -o /dev/stdout !{bam_abra} | samtools fixmate - - | sambamba sort -t !{half_threads} -m !{half_mem}G --tmpdir=sort_tmp -o "!{bam_tag}_sorted_fixmate.bam" /dev/stdin
    '''
}

