suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("-o", "--output",  
    dest="o", help="Name of the output table with defined ALGs 1-5.")

args <- parser$parse_args()

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
table_of_groups

dominant_alg_pairs <- strsplit(names(table_of_groups)[1:5], ' ')

for (ALG in 1:5){
    n1_name <- dominant_alg_pairs[[ALG]][1]
    n2_name <- dominant_alg_pairs[[ALG]][2]
    alg_rows <- which(syngraph[, 'n1'] == n1_name & syngraph[, 'n2'] == n2_name)

    chr_in_bibio <- names(sort(table(bibio_odb12[bibio_odb12[, 'busco_odb12'] %in% syngraph[alg_rows, 'busco'], 'chr']), T)[1])
    alg_label <- paste0('d', bibio_odb10_asn[chr_in_bibio, 'alg'])
    
    syngraph[alg_rows, 'ALG'] <- alg_label
}

### Comment out the next few lines if not to include the dot.
ALGs_pruned_with_dot <- read.table('data/ALG_assignments_pruned_m100.tsv', col.names = c('busco', 'ALG'))
ALG6_buscos <- ALGs_pruned_with_dot[ALGs_pruned_with_dot[, 'ALG'] == 'd6', ] # take only LG from the first good reconstruction of the ancestral state
ALG6_buscos <- ALG6_buscos[!(ALG6_buscos[, 'busco'] %in% syngraph[, 'busco']), ] # remove buscos assigned to other LGs

final_ALGs <- rbind(syngraph[!is.na(syngraph[, 'ALG']), c('busco', 'ALG')], ALG6_buscos)
final_ALGs <- final_ALGs[order(final_ALGs[, 2], final_ALGs[, 1]), ]

# 'data/syngraph.pruned2.100.ALGs.inferred.tsv'
write.table(final_ALGs, file = args$o, col.names = F, row.names = F, quote = F, sep = '\t')