suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("-o", "-output",  
    dest="o", help="Name of the output table with defined ALGs 1-5.")
parser$add_argument("-n", "-node", default = 'n1n2',  
    dest="n", help="Spedify node that defines the linkage groups (default, overlap of n1 and n2)")
parser$add_argument("-lgn", "-linkage_group_name", default = 'd',  
    dest="lgn", help="The prefix name for the specified lg")

args <- parser$parse_args()
# args$o <- 'tables/ALG_brachycera'
# args$n <- 'n21'
# args$lgn <- 'db'

# d - diptra
# db - brachycera
# ds - schizophora
# dm - culicidae (m for mosquitos)
# dc - chironomidae -> impossible with the current run

# read bibio BUSCOs;
# read old assignments
AGLs_odb10 <- read.table('data/Diptera_ALG_odb10.tsv', col.names = c('busco_odb10', 'ALG'))
bibio_odb10 <- read.table('data/Bibio_marci.buscos.odb10.tsv', col.names = c('busco_odb10', 'chr', 'from', 'to'))
odb10_labels_bibio <- merge(bibio_odb10, AGLs_odb10)
bibio_alt_assignment_table <- sort(table(paste(odb10_labels_bibio[ , 'chr'], odb10_labels_bibio[ , 'ALG'])), T)[1:5]

list_of_asn <- strsplit(names(bibio_alt_assignment_table), ' ')
bibio_odb10_asn <- data.frame(chr = sapply(list_of_asn, function(x){ x[1] } ), alg = sapply(list_of_asn, function(x){ x[2] } ))
row.names(bibio_odb10_asn) <- bibio_odb10_asn[, 'chr']
# this table contains chromosomal correspondences to ALGs in bibio. We used bibio because we know it has nice correspondences of the 5 chromosmes

bibio_odb12 <- read.table('data/busco_tables/Bibio_marci.syngraph.buscos.tsv', col.names = c('busco_odb12', 'chr', 'from', 'to'))

alg_groups <- list()

if (args$n == 'n1n2'){
    # read n1
    # read n2
    n1 <- read.table("data/syngraph/node_assignments/n1_asgn.tsv", col.names = c('busco', 'n1'))
    n2 <- read.table('data/syngraph/node_assignments/n2_asgn.tsv', col.names = c('busco', 'n2'))

    syngraph <- merge(n1, n2) # this will drop all those missing in one or the other
    syngraph[, 'ALG'] <- NA

    # double check both those are 5
    print(length(unique(n1[, 'n1'])))
    print(length(unique(n2[, 'n2'])))

    # sorted table shows 5 corresponding ALGs + several genes that moved around, those we can ignore
    table_of_groups <- sort(table(paste(syngraph[, 'n1'], syngraph[, 'n2'])), T) # 4 LGs, not 5
    print(table_of_groups)

    dominant_alg_pairs <- strsplit(names(table_of_groups)[1:5], ' ')


    for (ALG in 1:5){
        n1_name <- dominant_alg_pairs[[ALG]][1]
        n2_name <- dominant_alg_pairs[[ALG]][2]
        alg_rows <- which(syngraph[, 'n1'] == n1_name & syngraph[, 'n2'] == n2_name)

        alg_groups[[ALG]] <- syngraph[alg_rows, 'busco']
    }
} else {
    syngraph <- read.table(paste0("data/syngraph/node_assignments/", args$n ,"_asgn.tsv"), col.names = c('busco', 'node'))
    syngraph[, 'ALG'] <- NA

    table_of_groups <- sort(table(syngraph[, 'node']), T)
    print(table_of_groups)

    for ( ALG in 1:length(table_of_groups)){
        alg_name <- names(table_of_groups)[ALG]
        alg_rows <- which(syngraph[, 'node'] == alg_name)
        alg_groups[[ALG]] <- syngraph[alg_rows, 'busco']
    }

}

row.names(syngraph) <- syngraph[, 'busco']
node2lg <- rep('X', length(table_of_groups))
names(node2lg) <- names(table_of_groups)

if (args$lgn == 'd'){
    for (ALG in 1:5){
        chr_in_bibio <- names(sort(table(bibio_odb12[bibio_odb12[, 'busco_odb12'] %in% alg_groups[[ALG]], 'chr']), T)[1])
        alg_label <- paste0('d', bibio_odb10_asn[chr_in_bibio, 'alg'])
        
        syngraph[alg_groups[[ALG]], 'ALG'] <- alg_label
        print(paste(names(table_of_groups)[ALG], "assigned as", alg_label, "with", length(alg_groups[[ALG]]), "marker genes"))
    }
} else {
    for (ALG in 1:length(table_of_groups)){
        node_name <- names(table_of_groups)[ALG]
        ALG2Bibio <- sort(table(bibio_odb12[bibio_odb12[, 'busco_odb12'] %in% alg_groups[[ALG]], 'chr']), T)
        chr_in_bibio <- names(ALG2Bibio[1])
        alg_label <- paste0(args$lgn, bibio_odb10_asn[chr_in_bibio, 'alg'])
        print(paste(names(table_of_groups)[ALG], "assigned as", alg_label, "with", length(alg_groups[[ALG]]), "marker genes"))

        if ( alg_label %in% node2lg){
            node2lg[alg_label == node2lg] <- paste0(alg_label, 'a')
            node2lg[node_name] <- paste0(alg_label, 'b')
        } else {
            node2lg[node_name] <- alg_label
        }
    }
}

syngraph[, 'ALG'] <- node2lg[syngraph[, 'node']]
    # syngraph[alg_groups[[ALG]], 'ALG'] <- alg_label



### Comment out the next few lines if not to include the dot.
ALGs_pruned_with_dot <- read.table('data/ALG6_BUSCOs.tsv', col.names = c('busco', 'ALG'))
ALG6_buscos <- ALGs_pruned_with_dot[ALGs_pruned_with_dot[, 'ALG'] == 'd6', ] # take only LG from the first good reconstruction of the ancestral state
# sum((ALG6_buscos[, 'busco'] %in% syngraph[, 'busco']))
ALG6_buscos <- ALG6_buscos[!(ALG6_buscos[, 'busco'] %in% syngraph[, 'busco']), ] # remove buscos assigned to other LGs
if ( args$lgn != 'd'){
    ALG6_buscos$ALG <- paste0(args$lgn, 6)
}

final_ALGs <- rbind(syngraph[!is.na(syngraph[, 'ALG']), c('busco', 'ALG')], ALG6_buscos)
final_ALGs <- final_ALGs[order(final_ALGs[, 2], final_ALGs[, 1]), ]

# 'data/syngraph.pruned2.100.ALGs.inferred.tsv'
# args$o = 'tables/ALG_syngraph.pruned2.100.ALGs_brachycera.tsv'
write.table(final_ALGs, file = args$o, col.names = F, row.names = F, quote = F, sep = '\t')