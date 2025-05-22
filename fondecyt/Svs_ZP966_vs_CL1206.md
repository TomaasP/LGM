# Deteccion de SVs entre CL1206 y ZP966 a partir de genomas compeltos.

1. Llamar SVs con mum&Co 

`bash PathTo/mumandco_v3.8.sh -r .PathToRef/genome.fa -q PathToquery/Genome.fa -g 12500000 -o Query_Ref_Svs -ml 1000`  

Raw result:  
1206_966_Svs  Total SVs  =  184  
1206_966_Svs  Deletions  =  84  
1206_966_Svs  Insertions  =  83  
1206_966_Svs  Duplications  =  5  
1206_966_Svs  Contractions  =  3  
1206_966_Svs  Inversions  =  2  
1206_966_Svs  Translocations  =  7  

2. Filtrar por posicion y tamaño (conservar aquellas mayores a 1 kb remover translocaciones y eventos en el cromosoma XII y los eventos complejos)  

[chrXII fitlered](./chrXII_filtered.txt) -> 153 SVs  
[chrXII + > 1kb filtered](./chrXII_size_filtered.txt) -> 81  
[Trasloc](./map_966.assembly.final_to_1206.assembly.final.svg) -> 2 (ChrII_ChrVII, ChrXV_ChrII)  

3. Distribucion de SVs (en ambos genomas)  
[code:](./kariotypeR.r)  
[CL1206](./CL1206_Svs.pdf)    
[ZP966]()  

4. Densidad de genes en cada cromosoma
* Generar ventanas de 1 kb a partir del indexado de un genoma  

`bedtools makewindows -g 1206.assembly.final.fa.fai -w 1000 > 1206.assembly.final.1kbW.bed`  

* Preparar el archivo gff3 para que solo contenga las coordenadas de los genes.  

`awk '$3 == "gene" {print $1"\t"($4 - 1)"\t"$5}' CL1206.final.gff3 > CL1206.final.genes.bed`  

* Intersectar para sumar una cuerta columna con el numero de genes por ventana.  

`bedtools intersect -a 1206.assembly.final.1kbW.bed -b CL1206.final.genes.bed -c > 1206_geneDensity.txt`  
* Iterar el [](./kariotypeR.r)  

[**Resultado**](./chromosome.pdf)  

5. Genes de CL1206 con SVs respecto a ZP966  

* Procesar el gff para obtener un archivo con el nombre de cada cromosoma, el inicio, el final y el nombre de cada feature.  

`gawk '$3 == "gene" { start = $4 - 1 end = $5 chr = $1 match($9, /Name=([^;]+)/, name) print chr"\t"start"\t"end"\t"name[1]}' CL1206.final.gff3 > CL1206_genes_with_names.bed`  

* Intersect el archivo con la posicion de cada SVs y el de genes con nombres  

`bedtools intersect -a CL1206_genes_with_names.bed -b CL1206_Svs.bed -wa -wb > CL1206_vs_ZP966_genes_with_SVs.bed`  

[Resultado](./CL1206_vs_ZP966_genes_with_SVs.bed)  
Las primeras cuatro columnas corresponden a la info del gen, y las siguientes a la SVs detectada.  

6. Extraccion de genes y analisis de enriquecimiento Go  

`cut -f4 CL1206_vs_ZP966_genes_with_SVs.bed | sort | uniq > uniq_CL1206_vs_ZP966_genes_with_SVs.txt`  

[Resultado](./uniq_CL1206_vs_ZP966_genes_with_SVs.txt)  

* Analisis de GO en la plataforma gProfiler  

https://biit.cs.ut.ee/gplink/l/aj5NGZ1GXTm  

[Resultado](./gProfiler_scerevisiae_2025-05-22_10-54-47.png)  

[Tabla](./gProfiler_scerevisiae_22-5-2025_12-54-08%20p.m.__intersections.csv)  

Hay algunas asociaciones, pero tampoco un enriquecimiento tan notorio en algun fenotipo. Quizas, bien rebuscadmente, haya algo asociado a membranas.  






