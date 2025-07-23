log.info """\
         I N F E R   A L G   N F   P I P E L I N E    
         ===================================
         Input TSV : ${params.input_tsv}
         BUSCO DB : ${params.busco_db}
         BUSCO download path : ${params.busco_download_path}
         iqtree model : ${params.iqtree_model}
         iqtree outgroup : ${params.iqtree_outgroup}
         Taxon : ${params.taxon}
         outdir : ${params.outdir}
         """
         .stripIndent()

include { infer_alg_flow; to_busco; from_busco } from './infer_alg_flows.nf'

// full workflow
//workflow {
//         infer_alg_flow(
//                params.input_tsv, 
//                params.busco_db, 
//                params.busco_download_path, 
//                params.iqtree_model, 
//                params.iqtree_outgroup, 
//                params.syngraph_reference,
//                params.taxon
//                )
//}

// partitioned workflows
workflow {
       to_busco(
              params.input_tsv, 
              params.busco_db, 
              params.busco_download_path, 
              params.taxon
              )
}

//workflow{
//       from_busco(
//              params.from_busco_input_tsv, 
//              params.iqtree_model, 
//              params.iqtree_outgroup, 
//              params.syngraph_reference,
//              params.taxon
//              )
//}


// orthofinder version?

// synteny plotting?

// adding non ncbi genomes?