#! /usr/bin/env nextflow

// Copyright (C) 2025 IARC/WHO

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

nextflow.enable.dsl=2

// ---------------------------
// PARAMETERS
// ---------------------------

params.help = null
params.tumor_bam_folder = null
params.normal_bam_folder = null
params.bam_folder = null
params.bed = null
params.single = null
params.ref = null
params.abra_path = '/opt/conda/envs/abra-nf/share/abra2*/abra2.jar'
params.junctions = null
params.gtf = null
params.rna = null
params.ignore_bad_assembly = null
params.suffix_tumor  = "_T"
params.suffix_normal = "_N"
params.mem = 16
params.cpu = 4
params.output_folder = "abra_BAM"

log.info ""
log.info "--------------------------------------------------------"
log.info "  abra2-nf v4.0: Nextflow pipeline for ABRA2         "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

// ---------------------------
// HELP MESSAGE
// ---------------------------

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
    log.info '    --bed                FILE                    Bed file containing target intervals.'
    log.info '    --junctions                                  Flag to use STAR identified junctions.'
    log.info '    --gtf                FILE                    GTF file containing junction annotations.'
    log.info '    --rna                                        Flag to add RNA-specific recommended ABRA2 parameters.'
    log.info '    --mem                INTEGER                 RAM used (in GB, default: 16)'
    log.info '    --threads            INTEGER                 Number of threads (default: 4)'
    log.info '    --output_folder      FOLDER                  Output folder (default: abra_BAM).'
    log.info ''
    exit 0
}

 else {
      /* Software information */
  log.info "bam_folder = ${params.bam_folder}"
   log.info "ref          = ${params.ref}"
   log.info "cpu          = ${params.cpu}"
   log.info "mem          = ${params.mem}"
   log.info "output_folder= ${params.output_folder}"
   log.info "bed          = ${params.bed}"
   log.info "abra_path    = ${params.abra_path}"
   log.info "gtf          = ${params.gtf}"
   log.info "junctions    = ${params.junctions}"
   log.info "help=${params.help}"
 }


// ---------------------------
// PARAMETER CHECKS
// ---------------------------

assert params.ref, "please specify --ref option (--ref reference.fasta(.gz))"
assert params.abra_path, "please specify --abra_path option (--abra_path /path/to/abra.jar)"

if (params.bam_folder) {
    assert file(params.bam_folder).exists() : "BAM folder not found: ${params.bam_folder}"
} else {
    assert file(params.tumor_bam_folder).exists() : "Tumor BAM folder not found: ${params.tumor_bam_folder}"
    assert file(params.normal_bam_folder).exists() : "Normal BAM folder not found: ${params.normal_bam_folder}"
}

// ---------------------------
// FILE DEFINITIONS
// ---------------------------


fasta_ref = file(params.ref)
fasta_ref_fai = file( params.ref+'.fai' )
fasta_ref_sa  = file( params.ref+'.sa' )
fasta_ref_bwt = file( params.ref+'.bwt' )
fasta_ref_ann = file( params.ref+'.ann' )
fasta_ref_amb = file( params.ref+'.amb' )
fasta_ref_pac = file( params.ref+'.pac' )

bed = params.bed ? file(params.bed) : file('nothing')
gtf = params.gtf ? file(params.gtf) : file('nothing')

// ---------------------------
// PROCESSES
// ---------------------------

    process ABRA_SINGLE {
        tag { bam_tag }
        cpus params.cpu
        memory "${params.mem} GB"

        input:
		tuple val(bam_tag), path(bam), path(bai), path(junctions)
		path bed
		path fasta_ref
		path fasta_ref_fai
		path fasta_ref_sa
		path fasta_ref_bwt
		path fasta_ref_ann
		path fasta_ref_amb
		path fasta_ref_pac

        output:
		path("${bam_tag}_abra.ba*")

        script:
		java_mem = params.mem - 2
		abra_single   = params.single ? '--single --mapq 20' : ''
		abra_bed      = params.bed ? "--targets $bed" : ''
		abra_junctions= params.junctions ? "--junctions $junctions" : ''
		abra_gtf      = params.gtf ? "--gtf $gtf" : ''
		abra_rna      = params.rna ? "--sua --dist 500000" : ''
		abra_iba      = params.ignore_bad_assembly ? "--ignore-bad-assembly" : ''
    """
	    java -Xmx${java_mem}g -jar ${params.abra_path} \
			--in ${bam} --out ${bam_tag}_abra.bam \
			--ref ${fasta_ref} --tmpdir . \
			--threads ${params.cpu} --index ${abra_single} \
			${abra_bed} ${abra_junctions} ${abra_gtf} ${abra_rna} ${abra_iba} \
			> ${bam_tag}_abra.log 2>&1
    """
    }

	process ABRA_TN {
		tag { tumor_normal_tag }
		cpus params.cpu
		memory "${params.mem} GB"
		publishDir params.output_folder, mode: 'move'

		input:
		tuple val(tumor_normal_tag), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
		path bed
		path fasta_ref
		path fasta_ref_fai
		path fasta_ref_sa
		path fasta_ref_bwt
		path fasta_ref_ann
		path fasta_ref_amb
		path fasta_ref_pac

		output:
		path("${tumor_normal_tag}${params.suffix_normal}_abra.ba*")
		path("${tumor_normal_tag}${params.suffix_tumor}_abra.ba*")

		script:
		abra_single = params.single ? '--single --mapq 20' : ''
		abra_bed = params.bed ? "--targets $bed" : ''
    """
		java -Xmx${params.mem}g -jar ${params.abra_path} \
			--in ${normal_bam},${tumor_bam} \
			--out ${tumor_normal_tag}${params.suffix_normal}_abra.bam,${tumor_normal_tag}${params.suffix_tumor}_abra.bam \
			--ref ${fasta_ref} --threads ${params.cpu} --index ${abra_single} ${abra_bed} \
			> ${tumor_normal_tag}_abra.log 2>&1
    """
}

 // ---------------------------
// WORKFLOW
// ---------------------------

workflow {
    if (params.bam_folder) {
        log.info "Running single-sample ABRA2 realignment"
        bams = Channel.fromPath("${params.bam_folder}/*.bam")
            .map { bam -> tuple(bam.baseName, bam) }
        bais = Channel.fromPath("${params.bam_folder}/*.bam.bai")
            .map { bai -> tuple(bai.baseName.replace('.bai',''), bai) }

        if (params.junctions) {
            junctions = Channel.fromPath("${params.bam_folder}/*.SJ.out.tab")
                .map { j -> tuple(j.baseName.replace('.SJ.out.tab',''), j) }
            bam_bai = bams.join(bais).join(junctions)
        } else {
            bam_bai = bams.join(bais).map { tag, bam, bai -> tuple(tag, bam, bai, file("NO_JUNCTION_FILE")) }
        }

        ABRA_SINGLE(bam_bai, bed, fasta_ref, fasta_ref_fai, fasta_ref_sa, fasta_ref_bwt, fasta_ref_ann, fasta_ref_amb, fasta_ref_pac)

    } else {
        log.info "Running Tumor/Normal ABRA2 realignment"

        tumor_bams = Channel.fromPath("${params.tumor_bam_folder}/*${params.suffix_tumor}.bam")
            .map { bam -> tuple(bam.baseName.replace(params.suffix_tumor, ''), bam) }
        tumor_bais = Channel.fromPath("${params.tumor_bam_folder}/*${params.suffix_tumor}.bam.bai")
            .map { bai -> tuple(bai.baseName.replace(params.suffix_tumor, ''), bai) }
        tumor_bam_bai = tumor_bams.join(tumor_bais)

        normal_bams = Channel.fromPath("${params.normal_bam_folder}/*${params.suffix_normal}.bam")
            .map { bam -> tuple(bam.baseName.replace(params.suffix_normal, ''), bam) }
        normal_bais = Channel.fromPath("${params.normal_bam_folder}/*${params.suffix_normal}.bam.bai")
            .map { bai -> tuple(bai.baseName.replace(params.suffix_normal, ''), bai) }
        normal_bam_bai = normal_bams.join(normal_bais)

        tn_pairs = tumor_bam_bai.join(normal_bam_bai)
            .map { tag, tumor_bam, tumor_bai, normal_bam, normal_bai -> tuple(tag, tumor_bam, tumor_bai, normal_bam, normal_bai) }

        ABRA_TN(tn_pairs, bed, fasta_ref, fasta_ref_fai, fasta_ref_sa, fasta_ref_bwt, fasta_ref_ann, fasta_ref_amb, fasta_ref_pac)
    }
}
