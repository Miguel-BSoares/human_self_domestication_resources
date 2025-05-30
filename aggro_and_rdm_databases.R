##
library(biomaRt)
library(rGREAT)
library(BioMartGOGeneSets)

listDatasets(ensembl)
listEnsemblArchives()
listAttributes(mart)
listFilters(mart)
listOrganisms()
seqType
list
##aggro target genes
aggro_genes_id=c("ENSG00000189221","ENSG00000142319",
                "ENSG00000184845",
                "ENSG00000149295",
                "ENSG00000151577",
                "ENSG00000069696",
                "ENSG00000169676",
                "ENSG00000103546",
                "ENSG00000120907",
                "ENSG00000150594",
                "ENSG00000043591",
                "ENSG00000188778",
                "ENSG00000108576",
                "ENSG00000178394",
                "ENSG00000135312",
                "ENSG00000179546",
                "ENSG00000168830",
                "ENSG00000179097",
                "ENSG00000102468",
                "ENSG00000135914",
                "ENSG00000147246",
                "ENSG00000166736",
                "ENSG00000164270",
                "ENSG00000157219",
                "ENSG00000158748",
                "ENSG00000148680",
                "ENSG00000101200",
                "ENSG00000166148",
                "ENSG00000198049",
                "ENSG00000101405",
                "ENSG00000180914",
                "ENSG00000173261",
                "ENSG00000163873",
                "ENSG00000157764",
                "ENSG00000140945",
                "ENSG00000128573",
                "ENSG00000263001",
                "ENSG00000157103")

##random genes
random_genes_id=c("ENSG00000179195",
               "ENSG00000258691",
               "ENSG00000109519",
               "ENSG00000049283",
               "ENSG00000063761",
               "ENSG00000092295",
               "ENSG00000172482",
               "ENSG00000158092",
               "ENSG00000168591",
               "ENSG00000166262",
               "ENSG00000254536",
               "ENSG00000124279",
               "ENSG00000133121",
               "ENSG00000174564",
               "ENSG00000171204",
               "ENSG00000172421",
               "ENSG00000179403",
               "ENSG00000189114",
               "ENSG00000141750",
               "ENSG00000079332",
               "ENSG00000167699",
               "ENSG00000171103",
               "ENSG00000131828",
               "ENSG00000122674",
               "ENSG00000165233",
               "ENSG00000141698",
               "ENSG00000111860",
               "ENSG00000116095",
               "ENSG00000120662",
               "ENSG00000180815",
               "ENSG00000137809",
               "ENSG00000157349",
               "ENSG00000165801",
               "ENSG00000110104",
               "ENSG00000257057",
               "ENSG00000171533",
               "ENSG00000180354",
               "ENSG00000148339",
               "ENSG00000114942",
               "ENSG00000160111",
               "ENSG00000262165",
               "ENSG00000149133",
               "ENSG00000131503",
               "ENSG00000100385",
               "ENSG00000168938",
               "ENSG00000106624",
               "ENSG00000182095",
               "ENSG00000147127",
               "ENSG00000204103",
               "ENSG00000147378",
               "ENSG00000100373",
               "ENSG00000144659",
               "ENSG00000267680",
               "ENSG00000105287",
               "ENSG00000160951",
               "ENSG00000152256",
               "ENSG00000073111",
               "ENSG00000177096",
               "ENSG00000170275",
               "ENSG00000113595",
               "ENSG00000169016",
               "ENSG00000151151",
               "ENSG00000255398",
               "ENSG00000271425",
               "ENSG00000156140",
               "ENSG00000047634",
               "ENSG00000131389",
               "ENSG00000274209",
               "ENSG00000283149",
               "ENSG00000158764",
               "ENSG00000110497",
               "ENSG00000169403",
               "ENSG00000211452",
               "ENSG00000279961",
               "ENSG00000163539",
               "ENSG00000131650",
               "ENSG00000071991",
               "ENSG00000239704",
               "ENSG00000203965",
               "ENSG00000178573",
               "ENSG00000154978",
               "ENSG00000197479",
               "ENSG00000197991",
               "ENSG00000164919",
               "ENSG00000151929",
               "ENSG00000151704",
               "ENSG00000143768",
               "ENSG00000124493",
               "ENSG00000284431",
               "ENSG00000161642",
               "ENSG00000107937",
               "ENSG00000109685",
               "ENSG00000160654",
               "ENSG00000129460",
               "ENSG00000103067",
               "ENSG00000204599",
               "ENSG00000164440",
               "ENSG00000174501",
               "ENSG00000153815",
               "ENSG00000163947")

##old human genome database
old_human <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                    host="grch37.ensembl.org", 
                    path="/biomart/martservice" ,
                    dataset="hsapiens_gene_ensembl")

##new human genome database
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast") #changed mirror site, because otherwise it did not work

##############################################################################################
##############################################################################################
##old human genome aggro dataset
aggro_genes_old_human <- getSequence(id = aggro_genes_id,
                           type="ensembl_gene_id",
                           seqType="gene_exon_intron",
                           mart=old_human) 

aggro_genes_old_human_coords <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
                            mart = old_human, 
                            filters = "ensembl_gene_id", 
                            values = aggro_genes_id)

##new human genome aggro dataset
aggro_genes <- getSequence(id = aggro_genes_id,
            type="ensembl_gene_id",
            seqType="gene_exon_intron",
            mart=mart) 
exportFASTA(aggro_genes,"aggro_db") ## if needed

aggro_genes_coords <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
                        mart = mart, 
                        filters = "ensembl_gene_id", 
                        values = aggro_genes_id)

#################################
##old human genome random dataset
random_genes_old_human<- getSequence(id = random_genes_id,
                        type="ensembl_gene_id",
                        seqType="gene_exon_intron",
                        mart=old_human)

exportFASTA(random_genes_old_human,"random_old_db.fasta")

random_genes_old_human_coords <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
                             mart = old_human, 
                             filters = "ensembl_gene_id", 
                             values = random_genes_id)

##new human genome random dataset
random_genes<- getSequence(id = random_genes_id,
                           type="ensembl_gene_id",
                           seqType="gene_exon_intron",
                           mart=mart) 
exportFASTA(random_genes,"random_db.fasta")

random_genes_coords <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
                            mart = mart, 
                            filters = "ensembl_gene_id", 
                            values = random_genes_id)

### make table of old and new locations
locations_both_genomes <- cbind(aggro_genes_old_human_coords,aggro_genes_coords)

#########################################
###Other sequence types###


seqType="5utr" ## for 5-untranslated regions ##

random_genes_5utr<- getSequence(id = random_genes_id,
                                 type="ensembl_gene_id",
                                 seqType="5utr",
                                 mart=old_human)

random_genes_coords <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
                             mart = mart, 
                             filters = "ensembl_gene_id", 
                             values = random_genes_id)



seqType="coding_gene_flank" ## promoter regions ##
      upstream=100 ## associated with promoter regions

random_genes_flank<- getSequence(id = random_genes_id,
                           type="ensembl_gene_id",
                           seqType="coding_gene_flank",
                           upstream=100,
                           mart=old_human) 
