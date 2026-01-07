#! /usr/bin/env nextflow

// Copyright (C) 2026 IARC/WHO
// This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
// See the GNU General Public License for more details <http://www.gnu.org/licenses/>.

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

//Header for the IARC tools - logo generated using the following page : http://patorjk.com/software/taag  (ANSI logo generator)
def IARC_Header (){
     return  """
#################################################################################
# ██╗ █████╗ ██████╗  ██████╗██████╗ ██╗ ██████╗ ██╗███╗   ██╗███████╗ ██████╗  #
# ██║██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔═══██╗██║████╗  ██║██╔════╝██╔═══██╗ #
# ██║███████║██████╔╝██║     ██████╔╝██║██║   ██║██║██╔██╗ ██║█████╗  ██║   ██║ #
# ██║██╔══██║██╔══██╗██║     ██╔══██╗██║██║   ██║██║██║╚██╗██║██╔══╝  ██║   ██║ #
# ██║██║  ██║██║  ██║╚██████╗██████╔╝██║╚██████╔╝██║██║ ╚████║██║     ╚██████╔╝ #
# ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═════╝ ╚═╝ ╚═════╝ ╚═╝╚═╝  ╚═══╝╚═╝      ╚═════╝  #
# Nextflow pipelines for cancer genomics.########################################
"""
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

fasta_ref     = file(params.ref)
fasta_ref_fai = file("${params.ref}.fai")

bed = params.bed ? file(params.bed) : null
gtf = params.gtf ? file(params.gtf) : null

// ---------------------------
// PROCESSES
// ---------------------------

    process ABRA_SINGLE {
        tag { bam_tag }
        cpus params.cpu
        memory "${params.mem}GB"

        input:
		tuple val(bam_tag), path(bam), path(bai), path(junctions_file)
		path bed
		path fasta_ref
		path fasta_ref_fai

        output:
		path("${bam_tag}_abra.bam"), emit: bam_out
		path("${bam_tag}_abra.bai"), emit: bai_out
		path("${bam_tag}_abra.log"), emit: log_out

		publishDir params.output_folder, mode: 'move'

        script:
		def java_mem = params.mem - 2
		def threads_val = params.cpu ?: 1
		
		// Build ABRA2 options dynamically in Groovy
    	def abra_flags_list = []
// TO ADD: 		// if (params.single)                          abra_flags_list << "--single --mapq 20"
    	if (bed && bed.name != 'nothing')        abra_flags_list << "--targets ${bed}"
    	if (params.junctions && junctions_file && junctions_file.name != 'NO_JUNCTION_FILE')    abra_flags_list << "--junctions ${junctions_file}"
    	if (gtf && gtf.name != 'nothing')        abra_flags_list << "--gtf ${gtf}"
    	if (params.rna)                          abra_flags_list << "--sua --dist 500000"
    	if (params.ignore_bad_assembly)          abra_flags_list << "--ignore-bad-assembly"

    	def abra_flags = abra_flags_list.join(' ')

    """
	#!/bin/bash
    set -euo pipefail
	java -Xmx${java_mem}g -jar ${params.abra_path} \
        --in ${bam} \
        --out "${bam_tag}_abra.bam" \
        --ref ${fasta_ref} \
        --tmpdir . \
        --threads ${threads_val} \
        --index \
		--single --mapq 20 \
        ${abra_flags} > "${bam_tag}_abra.log" 2>&1
   """
    }

	process ABRA_TN {
		tag { sample_id }
		cpus params.cpu
		memory "${params.mem}GB"

		input:
		tuple val(sample_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
		path bed
		path fasta_ref
		path fasta_ref_fai

		output:
		path("${sample_id}${params.suffix_normal}_abra.bam"), emit: tumor_bam_out
		path("${sample_id}${params.suffix_normal}_abra.bai"), emit: tumor_bam_out
		path("${sample_id}${params.suffix_tumor}_abra.bam"), emit: normal_bam_out
		path("${sample_id}${params.suffix_tumor}_abra.bai"), emit: normal_bam_out
		path("${sample_id}_abra.log"), emit: log_out

		publishDir params.output_folder, mode: 'move'

		script:
		def java_mem = params.mem - 2
		def threads_val = params.cpu ?: 1
		
		// Build ABRA2 options dynamically in Groovy
    	def abra_flags_list = []
    	if (bed && bed.name != 'nothing')        abra_flags_list << "--targets ${bed}"
		// if (params.single)                          abra_flags_list << "--single --mapq 20"
    	// if (params.junctions && junctions_file && junctions_file.name != 'NO_JUNCTION_FILE')    abra_flags_list << "--junctions ${junctions_file}"
    	// if (gtf && gtf.name != 'nothing')        abra_flags_list << "--gtf ${gtf}"
    	// if (params.rna)                          abra_flags_list << "--sua --dist 500000"
    	// if (params.ignore_bad_assembly)          abra_flags_list << "--ignore-bad-assembly"

    	def abra_flags = abra_flags_list.join(' ')

    """
		java -Xmx${params.mem}g -jar ${params.abra_path} \
			--in ${normal_bam},${tumor_bam} \
			--out ${sample_id}${params.suffix_normal}_abra.bam,${sample_id}${params.suffix_tumor}_abra.bam \
			--ref ${fasta_ref} --threads ${threads_val} --index \
			 ${abra_flags} > "${sample_id}_abra.log" 2>&1
    """
}

 // ---------------------------
// WORKFLOW
// ---------------------------

workflow {

  		log.info IARC_Header()
// --------------------------------------------------
// INFO / HELP
// --------------------------------------------------

log.info ""
log.info "----------------------------------------------------------------------------------------------------------------"
log.info "  abra2-nf v4.0: Nextflow pipeline for ABRA2         "
log.info "----------------------------------------------------------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it under certain conditions; see LICENSE for details."
log.info "----------------------------------------------------------------------------------------------------------------"
log.info ""

if (params.help) {
    log.info ''
	log.info '-------------------------------------------------------------'
    log.info 'USAGE: '
    log.info 'nextflow run iarcbioinf/abra-nf --tumor_bam_folder tumor_BAM/ --normal_bam_folder normal_BAM/ --ref ref.fasta'
	log.info '-------------------------------------------------------------'
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

    if (params.bam_folder) {
        log.info "Running single-sample ABRA2 realignment"
        bams = Channel.fromPath("${params.bam_folder}/*.bam")
			.map { f -> tuple(f.baseName, f) }
			.ifEmpty { error "No BAM files found in ${params.bam_folder}" }
        bais = Channel.fromPath("${params.bam_folder}/*.bam.bai")
            .map { f -> tuple(f.baseName.replace('.bam',''), f) }
			.ifEmpty { error "No BAI files found in ${params.bam_folder}" }
        
		bam_bai = bams.join(bais) // emit tag, bam, bai

		// Attach junction file or NONE per BAM
        bam_bai = bam_bai.map { tag, bam, bai ->
            junction_file = params.junctions ? file("${params.bam_folder}/STAR.${tag}.SJ.out.tab") : file("NO_JUNCTION_FILE")
            if (!junction_file.exists()) junction_file = file("NO_JUNCTION_FILE")
            tuple(tag, bam, bai, junction_file)
        }

        //bam_bai.view { "BAM_BAI → $it" }

        ABRA_SINGLE(bam_bai, bed, fasta_ref, fasta_ref_fai)

    } else {
        log.info "Running Tumor/Normal ABRA2 realignment"

		tumor_bams = Channel.fromPath("${params.tumor_bam_folder}/*${params.suffix_tumor}.bam")
            .map { f -> tuple(f.baseName.replace(params.suffix_tumor, ''), f) }

        tumor_bais = Channel.fromPath("${params.tumor_bam_folder}/*${params.suffix_tumor}.bam.bai")
            .map { f -> tuple(f.baseName.replace("${params.suffix_tumor}.bam", ''), f) }

        normal_bams = Channel.fromPath("${params.normal_bam_folder}/*${params.suffix_normal}.bam")
            .map { f -> tuple(f.baseName.replace(params.suffix_normal, ''), f) }

        normal_bais = Channel.fromPath("${params.normal_bam_folder}/*${params.suffix_normal}.bam.bai")
            .map { f -> tuple(f.baseName.replace("${params.suffix_normal}.bam", ''), f) }

		tumor_bb = tumor_bams.join(tumor_bais)
        normal_bb = normal_bams.join(normal_bais)

    	tn_pairs = tumor_bb
            .join(normal_bb)
            .map { tag, t_bam, t_bai, n_bam, n_bai ->
                tuple(tag, t_bam, t_bai, n_bam, n_bai)
            }

        tn_pairs.view { "TN_PAIR → $it" }

        ABRA_TN(tn_pairs, bed, fasta_ref, fasta_ref_fai)
	    }
}
