#! /usr/bin/env nextflow

// Copyright (C) 2017 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help = null
params.tumor_bam_folder = null
params.normal_bam_folder = null
params.bam_folder = null
params.targets = null
params.single = null
params.ref = null
params.abra_path = null
params.junctions = null
params.gtf = null
params.rna = null
params.ignore_bad_assembly = null  

log.info ""
log.info "--------------------------------------------------------"
log.info "  abra2-nf v2.0: Nextflow pipeline for ABRA2         "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info '  USAGE              '
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
    log.info '    --ref                FILE (with index)       Reference fasta file indexed.'
    log.info '    --abra_path          FILE                    abra.jar explicit path.'
    log.info 'Optional arguments:'
    log.info '   When using Tumor/Normal pairs:'
    log.info '    --suffix_tumor       STRING                  Suffix identifying tumor bam (default: "_T").'
    log.info '    --suffix_normal      STRING                  Suffix identifying normal bam (default: "_N").'
    log.info '   In all cases:'
    log.info '    --single                                     Flag for single-end sequencing.'
    log.info '    --bed                FILE                    Bed file containing intervals.'
    log.info '    --junctions                                  Flag to use STAR identified junctions.'
    log.info '    --gtf                FILE                    GTF file containing junction annotations.'
    log.info '    --rna                                        Flag to add RNA-specific recommended ABRA2 parameters.'
    log.info '    --mem                INTEGER                 RAM used (in GB, default: 16)'
    log.info '    --threads            INTEGER                 Number of threads (default: 4)'
    log.info '    --output_folder      FOLDER                  Output folder (default: abra_BAM).'
    log.info ''
    exit 1
}else {
  /* Software information */
  log.info "bam_folder = ${params.bam_folder}"
  log.info "ref          = ${params.ref}"
  log.info "cpu          = ${params.cpu}"
  log.info "mem          = ${params.mem}"
  log.info "output_folder= ${params.output_folder}"
  log.info "targets      = ${params.targets}"
  log.info "abra_path    = ${params.abra_path}"
  log.info "gtf          = ${params.gtf}"
  log.info "junctions    = ${params.junctions}"
  log.info "help=${params.help}"
}


assert (params.ref != true) && (params.ref != null) : "please specify --ref option (--ref reference.fasta(.gz))"
if (params.bam_folder) {
    assert (params.bam_folder != true) && (params.bam_folder != null) : "please specify --bam_folder option (--bam_folder bamfolder)"
} else {
    assert (params.normal_bam_folder != true) && (params.normal_bam_folder != null) : "please specify --normal_bam_folder option (--normal_bam_folder bamfolder)"
    assert (params.tumor_bam_folder != true) && (params.tumor_bam_folder != null) : "please specify --tumor_bam_folder option (--tumor_bam_folder bamfolder)"
}

if (params.targets!=null) {
    assert (params.targets != true) : "please specify file when using --bed option (--bed regions.bed)"
    try { assert file(params.targets).exists() : "\n WARNING : input bed file not located in execution directory" } catch (AssertionError e) { println e.getMessage() }
}
targets = params.targets ? file(params.targets) : file('nothing')

if (params.gtf!=null) {
    assert (params.gtf != true) : "please specify file when using --gtf option (--gtf annotations.gtf)"
    try { assert file(params.gtf).exists() : "\n WARNING : input gtf file not located in execution directory" } catch (AssertionError e) { println e.getMessage() }
}
gtf = params.gtf ? file(params.gtf) : file('nothing')

assert (params.abra_path != true) && (params.abra_path != null) : "please specify --abra_path option (--abra_path /path/to/abra.jar)"

fasta_ref = file(params.ref)
fasta_ref_fai = file( params.ref+'.fai' )
fasta_ref_sa  = file( params.ref+'.sa' )
fasta_ref_bwt = file( params.ref+'.bwt' )
fasta_ref_ann = file( params.ref+'.ann' )
fasta_ref_amb = file( params.ref+'.amb' )
fasta_ref_pac = file( params.ref+'.pac' )

params.suffix_tumor  = "_T"
params.suffix_normal = "_N"
params.mem = 16
params.cpu = 4
params.output_folder = "abra_BAM"


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
	      //.println()

    // recovering of bai files
    bais = Channel.fromPath( params.bam_folder+'/*.bam.bai' )
              .ifEmpty { error "Cannot find any bai file in: ${params.bam_folder}" }
              .map {  path -> [ path.name.replace(".bam.bai",""), path ] }
              //.println()

    if(params.junctions){
	// recovering junctions files
    	junctions = Channel.fromPath( params.bam_folder+'/*.SJ.out.tab' )
              .ifEmpty { error "Cannot find any junction tab files in: ${params.bam_folder}" }
              .map {  path -> [ path.name.replace(".SJ.out.tab","").replace("STAR.",""), path ] }
	// building bam-bai pairs
        bam_bai = bams
              .join(bais)
	      .join(junctions)
    }else{
    // building bam-bai pairs
    bam_bai = bams
              .join(bais)
              //.map { bam, bai -> [ bam[1], bai[1] ] }
    }
    
    process abra {
        cpus params.cpu
        memory params.mem+'GB'

        tag { bam_tag }

        publishDir params.output_folder, mode: 'move'

        input:
        set bam_tag, file(bam), file(bai), file(junctions) from bam_bai
        file targets
        file fasta_ref
        file fasta_ref_fai
        file fasta_ref_sa
        file fasta_ref_bwt
        file fasta_ref_ann
        file fasta_ref_amb
        file fasta_ref_pac

        output:
        file("${bam_tag}_abra.ba*") into bam_out

        shell:
	java_mem = params.mem - 2
        abra_single = params.single ? '--single --mapq 20' : ''
        abra_targets = params.targets ? "--targets $targets" : ''
	abra_junctions = params.junctions ? "--junctions $junctions" : ''
	abra_gtf = params.gtf ? "--gtf $gtf" : ''
	abra_rna = params.rna ? "--sua --dist 500000" : ''
	abra_iba = params.ignore_bad_assembly ? "--ignore-bad-assembly ": ""
        '''
        java -Xmx!{java_mem}g -jar !{params.abra_path} --in !{bam_tag}.bam --out "!{bam_tag}_abra.bam" --ref !{fasta_ref} --tmpdir . --threads !{params.cpu} --index !{abra_single} !{abra_targets} !{abra_junctions} !{abra_gtf} !{abra_rna} > !{bam_tag}_abra.log 2>&1
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

        cpus params.cpu
        memory params.mem+'GB'

        tag { tumor_normal_tag }

        publishDir params.output_folder, mode: 'move'

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
        file("${tumor_normal_tag}${params.suffix_normal}_abra.ba*") into normal_output
        file("${tumor_normal_tag}${params.suffix_tumor}_abra.ba*") into tumor_output

        shell:
        tumor_normal_tag = tn[0].baseName.replace(params.suffix_tumor,"")
        abra_single = params.single ? '--single --mapq 20' : ''
        abra_bed = params.bed ? "--targets $bed" : ''
	      '''
        java -Xmx!{params.mem}g -jar !{params.abra_path} --in !{tumor_normal_tag}!{params.suffix_normal}.bam,!{tumor_normal_tag}!{params.suffix_tumor}.bam --out "!{tumor_normal_tag}!{params.suffix_normal}_abra.bam","!{tumor_normal_tag}!{params.suffix_tumor}_abra.bam" --ref !{fasta_ref} --threads !{params.cpu}  --index !{abra_single} !{abra_bed} > !{tumor_normal_tag}_abra.log 2>&1
	      '''
    }
}
