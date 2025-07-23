process get_taxon_info{
        scratch true

        input:
        tuple val(meta), val(taxid)

        output:
        tuple val(meta), path("${taxid}.taxonomy.tsv")

        script:
        """
        get_taxon_info.sh ${taxid}
        """
}

process collate_taxon_info{
        publishDir params.outdir, mode:'copy'
        scratch true

        input:
        val(meta)
        path(collected_taxonomy_tsvs, stageAs: "tsvs/*")

        output:
        path("${meta}.taxonomy_info.tsv")

        script:
        """
        cat tsvs/* > ${meta}.taxonomy_info.tsv
        """
}
process download_genomes{
        publishDir params.outdir, mode:'copy'
        scratch true

        input:
        val(meta)
        val(accessions_file)

        output:
        tuple val(meta), path("${meta}_genomes/*")

        script:
        """
        mkdir ${meta}_genomes
        mkdir ncbi
        datasets download genome accession --dehydrated --inputfile ${accessions_file}
        unzip ncbi_dataset.zip
        mv ncbi_dataset ncbi
        datasets rehydrate --directory ncbi
        mv ncbi/ncbi_dataset/data/*/* ${meta}_genomes
        """
}

process get_chromosome_names{
        publishDir params.outdir, mode:'copy'
        scratch true

        input:
        val(meta)
        val(accessions_file)

        output:
        tuple val(meta), path("${meta}_chromosome_files/*.chromosomes.txt"), emit: chromosomes
        tuple val(meta), path("${meta}_chromosome_files/*.x.txt"), emit: xs, optional: true
        tuple val(meta), path("${meta}_chromosome_files/*.y.txt"), emit: ys, optional: true

        script:
        """
        mkdir ${meta}_chromosome_files
        mkdir sequence_jsons

        parallel -j1 'datasets summary genome --report sequence accession {} > sequence_jsons/{}.json' :::: ${accessions_file}
        parallel -j1 'parse_json.py sequence_jsons/{}.json' :::: ${accessions_file}

        mv *.chromosomes.txt ${meta}_chromosome_files
        mv *.x.txt ${meta}_chromosome_files || echo "no x"
        mv *.y.txt ${meta}_chromosome_files || echo "no y"
        """
}

// need to add a retry clause for memory too high for default
process busco{
        publishDir params.outdir, mode:'copy'
        cpus 8
        memory '64G'
        scratch true

        input:
        tuple val(meta), path(genome)
        val(busco_db)
        path(busco_download_path)

        output:
        tuple val(meta), path("busco_results/${meta}/run_${busco_db}/busco_sequences"), emit: sequences
        tuple val(meta), path("busco_results/${meta}/run_${busco_db}/full_table.tsv"), emit: table
        tuple val(meta), path("busco_results/${meta}/short_summary.*.txt"), emit: summary

        script:
        """
        busco -i ${genome} -m genome -l ${busco_db} -c ${task.cpus} -f -o busco_results/${meta} --download_path ${busco_download_path} --metaeuk --offline --tar
        """
}

process busco2fasta{
        cpus 20
        memory '16G'
        scratch true

        input:
        val(meta)
        path(busco_results, stageAs: "busco_dirs/*")

        output:
        tuple val(meta), path("busco2fasta_results/*")

        script:
        """
        parallel -j${task.cpus} 'tar -C {//} -zxf {} ' ::: busco_dirs/*/*/busco_sequences/*.tar.gz
        busco2fasta.py -b busco_dirs -o busco2fasta_results -s protein -p 0.9
        """
}

process busco2fasta_full{
        cpus 20
        memory '16G'
        scratch true

        input:
        val(meta)
        path(busco_results, stageAs: "busco_dirs/*")

        output:
        tuple val(meta), path("busco2fasta_results/*")

        script:
        """
        parallel -j${task.cpus} 'tar -C {//} -zxf {} ' ::: busco_dirs/*/busco_sequences/*.tar.gz
        busco2fasta.py -b busco_dirs -o busco2fasta_results -s protein -p 0.9
        """
}

process mafft{
        publishDir params.outdir, mode:'copy'
        cpus 64
        memory '128G'
        queue 'week'

        input:
        tuple val(meta), path(fastas, stageAs: "*")

        output:
        tuple val(meta), path("mafft_alignments/*")

        script:
        """
        mkdir mafft_alignments
        mafft.sh ${task.cpus}
        """
}

process trimal{
        publishDir params.outdir, mode:'copy'
        input:
        tuple val(meta), path(fastas, stageAs: "*")

        output:
        tuple val(meta), path("trimmed_alignments/*")

        script:
        """
        mkdir trimmed_alignments
        trimal.sh
        """
}

// this needs to write new files
process trimal_array_clean{
        publishDir params.outdir, mode:'copy'
        cpus 8
        memory '3G'

        input:
        tuple val(meta), path(fastas, stageAs: "trimal/*")

        output:
        tuple val(meta), path("trimal_cleaned/*")

        script:
        """
        mkdir trimal_cleaned
        trimal_array_clean.py
        """
}

process catfasta2phyml{
        publishDir params.outdir, mode:'copy'
        cpus 16
        memory '10G'

        input:
        tuple val(meta), path(fastas, stageAs: "*")

        output:
        tuple val(meta), path("${meta}.supermatrix.phy"), emit: supermatrix
        tuple val(meta), path("${meta}.partitions.txt"), emit: partitions

        script:
        """
        catfasta2phyml.pl *.faa -c > ${meta}.supermatrix.phy 2> ${meta}.partitions.txt
        """
}

process iqtree{
        publishDir params.outdir, mode:'copy'
        cpus 160
        memory '320G'
        queue 'week'

        input:
        tuple val(meta), path(alignment)
        val(iqtree_model)
        val(iqtree_outgroup)

        output:
        tuple val(meta), path("*treefile")

        script:
        """
        iqtree -s ${alignment} -m "${iqtree_model}" -T ${task.cpus} -B 1000 -o ${iqtree_outgroup} 
        """
}

process prepare_busco_tables{
        publishDir "${params.outdir}/syngraph_busco_tables/", mode:'copy'

        input:
        tuple val(meta), path(busco_dir), path(chromosome_file)

        output:
        tuple val(meta), path("${meta}.syngraph.buscos.tsv")

        script:

        """
        cut -f 1,3,4,5 ${busco_dir}/run_*/full_table.tsv | sed '/^#/d' | awk '\$2 != ""' > ${meta}.busco.reformatted.tsv
        grep -f ${chromosome_file} ${meta}.busco.reformatted.tsv > ${meta}.syngraph.buscos.tsv 
        remove_blank_lines.sh ${meta}.syngraph.buscos.tsv
        """
}

process syngraph_build{
        publishDir params.outdir, mode:'copy'
        memory '24G'

        input:
        val(meta)
        path(input_dir, stageAs: "syngraph_busco_tables/*")

        output:
        tuple val(meta), path("${meta}.syngraph_build.pickle")

        script:
        """
        syngraph build -d syngraph_busco_tables -m -o ${meta}.syngraph_build
        """
}

process syngraph_infer{
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(meta), path(syngraph), path(tree)
        val(reference_taxon)

        output:
        tuple val(meta), path("${meta}.syngraph_infer*")

        script:
        """
        syngraph infer -g ${syngraph} -t ${tree} -m 89 -r 2 -a quick -s ${reference_taxon} -o ${meta}.syngraph_infer
        """
}

process syngraph_tabulate{
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(meta), path(syngraph_infer_files)

        output:
        tuple val(meta), path("${meta}.syngraph_tabulate*")

        script:
        """
        syngraph tabulate -g ${meta}.syngraph_infer.with_ancestors.pickle -o ${meta}.syngraph_tabulate
        """
}


process summarise_clusters{
        input:

        output:

        script:
        """
        echo NA
        """
}
