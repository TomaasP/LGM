# Paquetes
install.packages("RIdeogram")
install.packages("rsvg")
library(RIdeogram)
library(rsvg)

setwd("/Users/tomasp/Library/CloudStorage/Dropbox/Tomas_Tesis_Proyectos_LGM/Tesis/LongReadsProject/SVs-Related/SVs_tesis/08.Fondecyt_Fco")
# Datos de los cromosomas
chr_data <- read.delim("CL1206_kariotype.txt")

# Anotacion sitios de interés
site_data <- read.delim("CL1206_Svs.txt")

# Densidad de genes  
gene_data <- read.delim("1206_geneDensity.txt")

# Viz
ideogram(karyotype = chr_data, overlaid= gene_data, label = site_data, label_type= "marker", width=100,
         output = "kariotipo.svg"
         )
convertSVG("kariotipo.svg",device="pdf",
           #width=6,
           #height=8,
           #dpi=300
           )
