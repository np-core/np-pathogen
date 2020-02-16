#!/usr/bin/env nextflow

/*

Pipeline            np-sepsis
Version             v0.1

Description         Workflow for bacterial and viral read classification
                    from nanopore DNA / RNA for Queensland Genomics sepsis
                    projects.

Authors             Eike Steinig, Dan Rawlinson, Lachlan Coin

*/

// We should think carefully on how to design the configuration file, especially with Google Cloud in mind

log.info """
--------------------------------------
     NANOPORE SEPSIS WORKFLOW
--------------------------------------

outdir          : $params.outdir
observe         : $params.observe

dna             : $params.dna
barcodes        : $params.barcodes
quality         : $params.quality
length          : $params.length

- Controls -

batch_control    : null
negative_control : null

- Virus detection -

virus           : $params.virus
assembler       : $params.assembler
genome_size     : $params.genome_size
opts            : $params.opts
medaka_model    : $params.medaka_model
signature       : $params.signature

- RNA Nanopore -

- DNA Illumina -

illumina        : $params.illumina
assembler       : $params.illumina_assembler
opts            : $params.illumina_opts

- RNA Illumina -


--------------------------------------
"""


if (params.observe) {
     Channel
        .watchPath(params.observe)
        .map { file -> tuple(file.baseName, file) }
        .into { fastq_nanopore; nanostat_pre };

} else {
    Channel
        .fromPath(params.dna)
        .map { file -> tuple(file.baseName, file) }
        .into { fastq_nanopore; nanostat_pre };
}

process NanoStatPre {
    
    tag { id }
    label "ont"

    publishDir "$params.outdir/fastq", mode: "copy", pattern: "*.txt"

    input:
    set id, file(fq) from nanostat_pre

    output:
    file("${id}.prestats.txt")

    """
    NanoStat --threads $task.cpus --fastq $fq > ${id}.prestats.txt
    """

}

if (params.prefilter) {

    process NanoFilt {

        // TODO: make this optional

        tag { id }
        label "ont"

        publishDir "$params.outdir/fastq", mode: "copy", pattern: "*.txt"

        input:
        set id, file(fq) from fastq_nanopore

        output:
        set id, file("${id}.filtered.fq") into (nanostat_post, nanotax)

        """
        NanoFilt --quality $params.quality --length $params.length $fq > ${id}.filtered.fq
        """

    }

    process NanoStatPost {
        
        tag { id }
        label "ont"

        publishDir "$params.outdir/fastq", mode: "copy", pattern: "*.txt"

        input:
        set id, file(fq) from nanostat_post

        output:
        file("${id}.poststats.txt")

        """
        NanoStat --threads $task.cpus --fastq $fq > ${id}.poststats.txt
        """

    }

} else {
    fastq_nanopore.set { nanotax }
}


process Kraken2 {

    // TODO: Split processing into own rule (+ control)

    label "kraken2"
    publishDir "$params.outdir/kraken", mode: "copy", pattern: "*.tar.gz"

    input:
    set id, file(fq) from nanotax

    output:
    set id, file("${id}.out"), file("${id}.report") into nanopath_kraken_server
    file("${id}.tar.gz")

    """
    kraken2 --db $params.taxdb --threads $task.cpus --output ${id}.out --report ${id}.report --use-names $fq
    tar -cvzf ${id}.tar.gz ${id}.out ${id}.report
    """

}

process NanoPathKraken {

    label "ont"
    publishDir "$params.outdir", mode: "copy"
    
    input:
    set id, file(reads), file(report) from nanopath_kraken_server

    output:
    file("${id}.json")
    
    """
    nanopath server process-kraken -r $report -f $reads -o ${id}.json
    """

}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Section - Contamination flagging - Cadhla and I want to include a workflow that flags contamination
// and quantifies read abundances against a known spikein.

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Section - RNA reads - at the moment this is secondary, but because of viral identification actually kinda essential
// to include RNA analysis - we will test a DNA / RNA sequence procedure for a single flowcell with washing.

// For now this section will include a simple cDNA analysis workflow, but we need data first.

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Section - Viral - specific viral genome assembly and typing workflow, need to think carefully about RdRp tagging and
// further classfication.

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Metagenome assembly - not figured out yet, especially for viral stuff, but a simple one with polishing here:

// process Assembly {
    
//   tag { "$id:$params.assembler" }
//   label "assembly"

//   publishDir "$params.outdir/assembly", mode: "copy"

//   input:
//   set id, file(fq) from assembly_fastq

//   output:
//   file("${id}.assembly_info.txt")
//   set id, file("${id}.assembly.fasta"), file(fq) into (racon_assembly, pilon_assembly)
//   set id, file("${id}.assembly_graph.gfa") into bandage_assembly
  
//   script:

//   if ( params.assembler == 'flye' )

//     """
//     flye --nano-raw $fq --genome-size $params.genome_size $params.opts -t $task.cpus -o assembly
//     mv assembly/assembly_info.txt ${id}.assembly_info.txt
//     mv assembly/assembly.fasta ${id}.assembly.fasta
//     mv assembly/assembly_graph.gfa ${id}.assembly_graph.gfa
//     """

// }

// process Racon {

//     tag { id }
//     label "racon"

//     input:
//     set id, file(assembly), file(fastq) from racon_assembly

//     output:
//     set id, file("${id}.racon.fasta"), file(fastq) into medaka_racon

//     script:
//     """
//     minimap2 -x map-ont -t $task.cpus $assembly $fastq > assembly_1.paf
//     racon -m 8 -x -6 -g -8 -w 500 -t $task.cpus $fastq assembly_1.paf $assembly > assembly_consensus_1.fasta
//     minimap2 -x map-ont -t $task.cpus assembly_consensus_1.fasta $fastq > assembly_2.paf
//     racon -m 8 -x -6 -g -8 -w 500 -t $task.cpus $fastq assembly_2.paf assembly_consensus_1.fasta > assembly_consensus_2.fasta
//     minimap2 -x map-ont -t $task.cpus assembly_consensus_2.fasta $fastq > assembly_3.paf
//     racon -m 8 -x -6 -g -8 -w 500 -t $task.cpus $fastq assembly_3.paf assembly_consensus_2.fasta > assembly_consensus_3.fasta
//     minimap2 -x map-ont -t $task.cpus assembly_consensus_3.fasta $fastq > assembly_4.paf
//     racon -m 8 -x -6 -g -8 -w 500 -t $task.cpus $fastq assembly_4.paf assembly_consensus_3.fasta > ${id}.racon.fasta
//     """
// }

// process Medaka {

//     tag { id }
//     label "medaka"

//     publishDir "$params.outdir/medaka", mode: "copy"

//     input:
//     set id, file(racon_assembly), file(fastq) from medaka_racon

//     output:
//     set id, file("${id}.consensus.fasta") into medaka_consensus

//     """
//     medaka_consensus -i $fastq -d $racon_assembly -o racon_medaka -t $task.cpus -m $params.medaka_model
//     mv racon_medaka/consensus.fasta ./${id}.consensus.fasta
//     """

// }

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Short read section to support assembly by polishing, we don't expect Illumina data for now
// but could become available post-hoc (for research), so some support should be included

if (params.illumina) {

    // Pair assemblies with parsed Illumina files by ID:


    illumina_reads = Channel.fromFilePairs(params.illumina, flat: true)

    process Trimmomatic {

        tag { rid }
        label "trimmomatic"

        input:
        set rid, file(forward), file(reverse) from illumina_reads

        output:
        set rid, file("${rid}_1P.fq.gz"), file("${rid}_2P.fq.gz") into (pilon_illumina, shovill_illumina, pilon_assembly_illumina)

        """
        trimmomatic PE $forward $reverse \
        -threads $task.cpus -phred33 -baseout ${rid}.fq.gz \
        ILLUMINACLIP:$baseDir/resources/trimmomatic/all_adapters.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

    }
    
    medaka_pilon = medaka_consensus.cross(pilon_illumina)
                    .map { crossed ->
                        if (crossed[0][0] == crossed[1][0]){
                            tuple( crossed[0][0], crossed[0][1],crossed[1][1], crossed[1][2] )
                        } else {
                            null
                        }
                    }
                    .filter { it != null }
    
    process PilonCorrection {

        tag { aid }
        label "pilon"

        publishDir "$params.outdir/pilon_medaka", mode: "copy"

        input:
        set aid, file(consensus_assembly), file(forward), file(reverse) from medaka_pilon

        output:
        file("${aid}.medaka.pilon.fasta")


        """
        minimap2 -ax sr $consensus_assembly $forward $reverse > aln.sam && \
            samtools view -S -b aln.sam > aln.bam && samtools sort aln.bam -o alignment1.bam && \
            samtools index alignment1.bam

        pilon --genome $consensus_assembly --frags alignment1.bam --outdir correction1 --changes
        
        minimap2 -ax sr correction1/pilon.fasta $forward $reverse > aln.sam && \
            samtools view -S -b aln.sam > aln.bam && samtools sort aln.bam -o alignment2.bam && \
            samtools index alignment2.bam

        pilon --genome correction1/pilon.fasta  --frags alignment2.bam --outdir correction2 --changes && \
        mv correction2/pilon.fasta ${aid}.medaka.pilon.fasta
        """

    }

    process Shovill {

        tag { rid }
        label "shovill"

        publishDir "$params.outdir/$params.illumina_assembler", mode: "copy"

        input:
        set rid, file(forward), file(reverse) from shovill_illumina

        output:
        file("${rid}.fasta")

        """
        shovill --R1 ${forward} --R2 ${reverse} --cpus $task.cpus --ram $task.memory \
        --depth 100 --assembler $params.illumina_assembler --outdir $rid --force
        mv ${rid}/contigs.fa ${rid}.fasta
        """

    }   
    
    assembly_pilon = pilon_assembly.cross(pilon_assembly_illumina)
                    .map { crossed ->
                        if (crossed[0][0] == crossed[1][0]){
                            tuple( crossed[0][0], crossed[0][1], crossed[1][1], crossed[1][2] )
                        } else {
                            null
                        }
                    }
                    .filter { it != null }
    
    process AssemblyPilon {

        tag { aid }
        label "pilon"

        publishDir "$params.outdir/pilon_assembly", mode: "copy"

        input:
        set aid, file(assembly), file(forward), file(reverse) from assembly_pilon

        output:
        file("${aid}.assembly.pilon.fasta")


        """
        minimap2 -ax sr $assembly $forward $reverse > aln.sam && \
            samtools view -S -b aln.sam > aln.bam && samtools sort aln.bam -o alignment1.bam && \
            samtools index alignment1.bam

        pilon --genome $assembly --frags alignment1.bam --outdir correction1 --changes
        
        minimap2 -ax sr correction1/pilon.fasta $forward $reverse > aln.sam && \
            samtools view -S -b aln.sam > aln.bam && samtools sort aln.bam -o alignment2.bam && \
            samtools index alignment2.bam

        pilon --genome correction1/pilon.fasta  --frags alignment2.bam --outdir correction2 --changes && \
        mv correction2/pilon.fasta ${aid}.assembly.pilon.fasta
        """

    }

}