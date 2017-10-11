#!/usr/bin/env nextflow

// requires (in path):
// java

params.help = null
params.tumor_bam_folder = null
params.normal_bam_folder = null
params.bam_folder = null
params.bed = null
params.single = null
params.ref = null
params.abra_path = null

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info '  NEXTFLOW pipeline for ABRA2              '
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run iarcbioinf/abra-nf --tumor_bam_folder tumor_BAM/ --normal_bam_folder normal_BAM/ --ref ref.fasta'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '   When using Tumor/Normal pairs:'
    log.info '    --tumor_bam_folder   FOLDER                  Folder containing tumor BAM files.'
    log.info '    --normal_bam_folder  FOLDER                  Folder containing matched normal BAM files.'
    log.info '   In other cases:'
    log.info '    --bam_folder         FOLDER                  Folder containing BAM files.'
    log.info '   In all cases:'
    log.info '    --ref                FILE (with index)       Reference fasta file indexed by bwa.'
    log.info '    --abra_path          FILE                    abra.jar explicit path.'
    log.info 'Optional arguments:'
    log.info '   When using Tumor/Normal pairs:'
    log.info '    --suffix_tumor       STRING                  Suffix identifying tumor bam (default: "_T").'
    log.info '    --suffix_normal      STRING                  Suffix identifying normal bam (default: "_N").'
    log.info '   In all cases:'
    log.info '    --single                                     Flag for single-end sequencing.'
    log.info '    --bed                FILE                    Bed file containing intervals.'
    log.info '    --mem                INTEGER                 RAM used (in GB, default: 16)'
    log.info '    --threads            INTEGER                 Number of threads (default: 4)'
    log.info '    --out_folder         FOLDER                  Output folder (default: abra_BAM).'
    log.info ''
    exit 1
}

assert (params.ref != true) && (params.ref != null) : "please specify --ref option (--ref reference.fasta(.gz))"

if (params.bam_folder) {
    assert (params.bam_folder != true) && (params.bam_folder != null) : "please specify --bam_folder option (--bam_folder bamfolder)"
} else {
    assert (params.normal_bam_folder != true) && (params.normal_bam_folder != null) : "please specify --normal_bam_folder option (--normal_bam_folder bamfolder)"
    assert (params.tumor_bam_folder != true) && (params.tumor_bam_folder != null) : "please specify --tumor_bam_folder option (--tumor_bam_folder bamfolder)"
}

if (params.bed!=null) {
    assert (params.bed != true) : "please specify file when using --bed option (--bed regions.bed)"
    try { assert file(params.bed).exists() : "\n WARNING : input bed file not located in execution directory" } catch (AssertionError e) { println e.getMessage() }
}
bed = params.bed ? file(params.bed) : file('nothing')

assert (params.abra_path != true) && (params.abra_path != null) : "please specify --abra_path option (--abra_path /path/to/abra.jar)"

fasta_ref = file(params.ref)
fasta_ref_fai = file( params.ref+'.fai' )
fasta_ref_gzi = file( params.ref+'.gzi' )
fasta_ref_sa = file( params.ref+'.sa' )
fasta_ref_bwt = file( params.ref+'.bwt' )
fasta_ref_ann = file( params.ref+'.ann' )
fasta_ref_amb = file( params.ref+'.amb' )
fasta_ref_pac = file( params.ref+'.pac' )

params.suffix_tumor = "_T"
params.suffix_normal = "_N"
params.mem = 16
params.threads = 4
params.out_folder = "abra_BAM"

try { assert fasta_ref.exists() : "\n WARNING : fasta reference not located in execution directory. Make sure reference index is in the same folder as fasta reference" } catch (AssertionError e) { println e.getMessage() }
if (fasta_ref.exists()) {assert fasta_ref_fai.exists() : "input fasta reference does not seem to have a .fai index (use samtools faidx)"}

if (fasta_ref.exists() && params.ref.tokenize('.')[-1] == 'gz') {assert fasta_ref_gzi.exists() : "input gz fasta reference does not seem to have a .gzi index (use samtools faidx)"}

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

        publishDir params.out_folder, mode: 'move'

        input:
        file bam_bai
        file bed
        file fasta_ref
        file fasta_ref_fai
        file fasta_ref_gzi
        file fasta_ref_sa
        file fasta_ref_bwt
        file fasta_ref_ann
        file fasta_ref_amb
        file fasta_ref_pac

        output:
        file("${bam_tag}_abra.bam*") into bam_out

        shell:
        bam_tag = bam_bai[0].baseName
        abra_single = params.single ? '--single --mapq 20' : ''
        abra_bed = params.bed ? "--bed $bed" : ''
        '''
        java -Xmx!{params.mem}g -jar !{params.abra_path} --in !{bam_tag}.bam --out "!{bam_tag}_abra.bam" --ref !{fasta_ref} --tmpdir . --threads !{params.threads} --index !{abra_single} !{abra_bed} > !{bam_tag}_abra.log 2>&1
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

        publishDir params.out_folder, mode: 'move'

        input:
        file tn from tn_bambai
        file bed
        file fasta_ref
        file fasta_ref_fai
        file fasta_ref_gzi
        file fasta_ref_sa
        file fasta_ref_bwt
        file fasta_ref_ann
        file fasta_ref_amb
        file fasta_ref_pac

        output:
        file("${tumor_normal_tag}${params.suffix_normal}_abra.bam*") into normal_output
        file("${tumor_normal_tag}${params.suffix_tumor}_abra.bam*") into tumor_output

        shell:
        tumor_normal_tag = tn[0].baseName.replace(params.suffix_tumor,"")
        abra_single = params.single ? '--single --mapq 20' : ''
        abra_bed = params.bed ? "--bed $bed" : ''
	      '''
        java -Xmx!{params.mem}g -jar !{params.abra_path} --in !{tumor_normal_tag}!{params.suffix_normal}.bam,!{tumor_normal_tag}!{params.suffix_tumor}.bam --out "!{tumor_normal_tag}!{params.suffix_normal}_abra.bam","!{tumor_normal_tag}!{params.suffix_tumor}_abra.bam" --ref !{fasta_ref} --threads !{params.threads}  --index !{abra_single} !{abra_bed} > !{tumor_normal_tag}_abra.log 2>&1
	      '''
    }
}
