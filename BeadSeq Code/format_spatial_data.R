setwd("/Users/rstickel/Documents/mac_lab/slideseq/data/seq_data/puck_seq/bijective_maps")
bijective = read.table('180430_6.csv', header= T,sep = ',', row.names = 1)


bc_locations <- read.table('180413_7_barcodelocations.csv', header = T, sep = ',', row.names = 1)
bc_locations_bijective = bc_locations[row.names(bijective),]


