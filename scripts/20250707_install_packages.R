conditional_install <- function(package){
    if(!(package %in% installed.packages()[,"Package"])){
        install.packages(package, repos = "https://mirrors.nic.cz/R") # this mirror can be adjusted to any...
    } else {
        print(paste(package, 'is installed already.'))
    }
}

## CRAN packages
conditional_install('ape')
conditional_install('phytools')
conditional_install('rphylopic')

## for bioconductor packages
if(!("ggtree" %in% installed.packages()[,"Package"])){
    conditional_install('BiocManager')
    BiocManager::install("ggtree")
} else {
    print('ggtree is installed already.')
}