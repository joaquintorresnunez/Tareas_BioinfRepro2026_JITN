# Tarea Sesión 2.3

Bioinformática y análisis genómicos reproducibles

Joaquín Torres Nuñez

Profesor: Ricardo Verdugo

14/10/2026

---


## Análisis de control de calidad

1. ¿Cómo se llaman los archivos que contienen las tasas de datos perdidos por SNP y por muestra?

PLINK --miss

2. ¿Cuántas variantes se eliminaron por tener una tasa de datos perdidos mayor a 0.2?

Se eliminaron 4,680 variantes (SNPs) por tener una tasa de datos perdidos mayor al 20% (--geno 0.2).

Desglose de los resultados de tu ejecución:

    Variantes iniciales: 813,366

    Variantes eliminadas: 4,680

    Variantes finales retenidas: 808,686

    Muestras conservadas: 48 (no se aplicó filtrado de individuos --mind en este paso)

    Tasa global de genotipado: 98.98%

3. ¿Cuántos individuos tenían una tasa de datos perdidos mayor a 0.02?

Se eliminaron 0 personas (muestras) al aplicar el filtro de tasa de datos perdidos del 2% (--mind 0.02). Todas las muestras pasaron el control de calidad.

Resumen del resultado de la ejecución:

    Muestras iniciales: 48

    Muestras eliminadas: 0

    Tasa total de genotipado: 100% (Total genotyping rate is exactly 1)

    Variantes retenidas: 574,624 de 574,624

    Muestras finales retenidas: 48

4. Basados en los histogramas y en sus cálculos, ¿qué valores umbrales de datos perdidos para muestras y SNPs sugeriría?
0.02 (2%). Esto demuestra que la inmensa mayoría del array o panel tiene excelente cobertura. Un filtro más estricto de 0.02 removerá las variantes ruidosas o con fallos de hibridación/llamado de alelos sin perder número de individuos,  lo cual es ideal para análisis poblacionales (PCA, estructura, LD pruning).

## Revisión discrepancia de sexo

1. ¿Cuántos individuos fueron eliminados por discrepancia de sexo?
Corroboración  de sexos: Se detectaron 3 discrpancias, se eliminaron 3 individuos por discrpancias de sexo:
ARI022 ARI022
CDSJ176 CDSJ176
CDSJ118 CDSJ118

2.  ¿Qué riesgo(s) se corre(n) si no se eliminaran?
No eliminar las muestras con discrepancia de sexo compromete severamente la validez de todo el estudio genómico. Los principales riesgos serían la geneación falsos positivos en estudios de asociación, PCA, parentesco, etc.

   
## Filtrado de SNPs

- Genere un archvio con SNP autosomales
```
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' chilean_all48_hg19_6.bim > snp_1_22.txt
plink --bfile chilean_all48_hg19_6 --extract snp_1_22.txt --make-bed --out chilean_all48_hg19_7
```
```
Samples: 45 total (26 males, 19 females; all 45 classified as founders).
Variants: Reduced from 574,624 to 557,922 variants using --extract (16,702 variants removed).
Output generated: chilean_all48_hg19_7.bed, chilean_all48_hg19_7.bim, chilean_all48_hg19_7.fam.
```

Elimine los SNPs monomórficos (sin variabilidad). En el set hay, en este punto, 45 individuos, por lo tanto, la menor frecuencia alélica posible para un polimorfismo es 1/(2x45)=0,011. Utilizaremos este valor como mínima MAF. *** Nota: *** En el tutorial de GWAS además se borran los SNPs con frecuencia alélica menor a 0,05. Sin embargo eso no es necesario para estudios de estructura.

```
plink --bfile chilean_all48_hg19_7 --maf 0.011 --make-bed --out chilean_all48_hg19_8

98544 variants removed due to minor allele threshold(s)
```
Tarea

1. ¿Cuál es el nombre del primer conjunto de datos que solo contiene SNPs en autosomas? 

Se filtra los cromosomas autosomicos (1-22) y se redirige el aricho de entrada BIM a un archivo de texto plano snp_1_22.txt

2. ¿Cuántos SNPs se encontraban en cromosomas sexuales?   

574624-557922= al excluir cromosomas sexuales 16702 SNPs

3. ¿Como calcularía el número de cromosomas que porta cada uno de los alelos para cada SNP?

Se podria calcular el numero de cromosomas que portan los allelos de frecuencia mninma multiplincado su frecuencia por el numero total de cromocomas (90)
El numero de cromosomas que porta el alelo de mayor frecuencia seria el total de cromosomas- el numero ttoal de cromosomas que porta MAF.

## Borrar SNPs por filtro de HWE

Los SNPs cuyos genotipos tengan frecuencias que se desvíen demasiado desde los valores esperados en Equilibrio de Hardy-Weinberg (HWE) pueden tener errores de genotipificación y deben ser eliminados.
```
plink --bfile chilean_all48_hg19_8 --hardy

459378 variants loaded from .bim file.
```
```
awk '{ if ($9 <0.000001) print $0 }' plink.hwe > plinkzoomhwe.hwe
```

Filtre los SNPs que no cumplen con el filtro MAX>0,000001.

```
plink --bfile chilean_all48_hg19_8 --hwe 1e-6 --hwe-all --make-bed --out chilean_all48_hg19_9

Calculating allele frequencies... done.
--hwe: 1281 variants removed due to Hardy-Weinberg exact test.
458097 variants and 45 people pass filters and QC.
Among remaining phenotypes, 0 are cases and 45 are controls.
```

1. ¿Cuál es el nombre del archivo con los resultados de la prueba de HWE?

El nombre del archivo es plink.hwe
2. ¿Basándose en la distribución de los valores de p, le parece el umbral usado razonable o propondría otro valor?

Por qué es adecuado: En un escaneo genómico, la gran mayoría de los SNPs están en equilibrio. El umbral utilizado es aceptable puesto que engloba gran parte de posibles fallas tecncias que ocurrieron al momento de la secuanción detectadas por filtro HWE.


## Eliminar parentescos desconocidos

Las pruebas de parentesco asumen que los SNP no están correlacionados, vale decir que no tienen un fuerte ligamiento. Para generar una lista de SNPs sin correlaciones muy altas, exploremos regiones de inversión (inversion.txt [regiones con alto LD]) y eliminaremos SNPs con el argumento command --indep-pairwise.

```
plink --bfile chilean_all48_hg19_9 --exclude $T/inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
```
Asumiendo una muestra poblacional aleatoria, vamos a excluir todos los individuos que tengan un valor pihat ≥ 0,2.
```
plink --bfile chilean_all48_hg19_9 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2

```
```
head pihat_min0.2.genome | cut -c 1-77


   FID1     IID1     FID2     IID2 RT    EZ      Z0      Z1      Z2  PI_HAT
  CDSJ177  CDSJ177   ARI001   ARI001 UN    NA  0.4733  0.5267  0.0000  0.2634
  CDSJ021  CDSJ021   ARI001   ARI001 UN    NA  0.5359  0.4641  0.0000  0.2320
   ARI006   ARI006   ARI001   ARI001 UN    NA  0.4314  0.5686  0.0000  0.2843
   ARI021   ARI021   ARI018   ARI018 UN    NA  0.6915  0.1554  0.1530  0.2307
   ARI021   ARI021   ARI015   ARI015 UN    NA  0.6506  0.2759  0.0736  0.2115
   ARI021   ARI021   ARI001   ARI001 UN    NA  0.4205  0.5795  0.0000  0.2897
   ARI021   ARI021   ARI019   ARI019 UN    NA  0.6748  0.2174  0.1078  0.2165
  CDSJ174  CDSJ174   ARI001   ARI001 UN    NA  0.5196  0.4804  0.0000  0.2402
  CDSJ175  CDSJ175   ARI001   ARI001 UN    NA  0.4559  0.5441  0.0000  0.2720

Rscript --no-save $T/Relatedness.R
```

```
Usemos R para investigar cuáles son los individuos que están generando más problemas.
> rel <- read.table("pihat_min0.2.genome", header=T)
> table(rel$IID1)

 ARI001  ARI002  ARI003  ARI004  ARI006  ARI008  ARI013  ARI015  ARI017  ARI018 
     11       1       1       1       1       1       1       1       1       2 
 ARI021  ARI023 CDSJ020 CDSJ021 CDSJ046 CDSJ099 CDSJ122 CDSJ157 CDSJ174 CDSJ175 
      4       1       1       1       1       1       1       1       1       1 
CDSJ177 CDSJ262 CDSJ283 CDSJ284 CDSJ297 CDSJ299 CDSJ347 CDSJ361 CDSJ372 CDSJ400 
      1       1       1       1       1       1       1       1       1       1 
CDSJ403 CDSJ469 CDSJ471 CDSJ472 
      1       1       1       1 

```


Claramente el individuo ARI001 genera la mayoría de los aparentes parentescos. Revisemos qué parentescos quedan si eliminamos este 
individuo.

```
> subset(rel, !IID1 %in% "ARI001" & !IID2 %in% "ARI001")
     FID1   IID1   FID2   IID2 RT EZ     Z0     Z1     Z2 PI_HAT PHE      DST
4  ARI021 ARI021 ARI018 ARI018 UN NA 0.6915 0.1554 0.1530 0.2307  -1 0.848091
5  ARI021 ARI021 ARI015 ARI015 UN NA 0.6506 0.2759 0.0736 0.2115  -1 0.841702
7  ARI021 ARI021 ARI019 ARI019 UN NA 0.6748 0.2174 0.1078 0.2165  -1 0.843969
16 ARI018 ARI018 ARI014 ARI014 UN NA 0.4936 0.4438 0.0626 0.2845  -1 0.852064
48 ARI008 ARI008 ARI019 ARI019 UN NA 0.6189 0.3611 0.0200 0.2005  -1 0.837686
   PPC  RATIO
4    1 2.6656
5    1 3.1995
7    1 2.8519
16   1 4.0182
48   1 3.3516
```
```
awk '$2=="ARI001" || $2=="ARI021" || $2=="ARI018"' chilean_all48_hg19_9.fam > to_romeve_by_relatedness.txt
```

Tambien se elimina ARI21 Y 18. 

```
plink -bfile chilean_all48_hg19_9 -remove to_romeve_by_relatedness.txt -make-bed --out chilean_all48_hg19_10

Calculating allele frequencies... done.
458097 variants and 42 people pass filters and QC.
Among remaining phenotypes, 0 are cases and 42 are controls.
--make-bed to chilean_all48_hg19_10.bed + chilean_all48_hg19_10.bim +
chilean_all48_hg19_10.fam ... done.
```

Eliminacion de 3 individuos mas, quedan 42.

1. ¿Cuántos SNPs en aparente equilibrio de ligamiento se encontraron?

wc -l indepSNP.prune.in
346968 of 450182 variants removed.
103214 SNPS en aparente equilibrio

2. ¿Cuántos SNPs se eliminaron por estar en regiones de inversiones conocidas?

Se eliminaron 7,915 SNPs por estar en las regiones de alto LD o inversiones conocidas especificadas en tu comando.

3. ¿Cuántos individuos quedaron luego del filtro de parentesco?

Se elminaron 3 nuevos individuos por filtro de parentesco, dejando 42 individuos en el estudio

4. ¿Cuál fue el mayor coeficiente de parentesco efectivamente aceptado?

PI_HAT (o π̂) representa la proporción estimada de identidad compartida por descendencia (IBD) entre un par de individuos, utilizada para medir el parentesco genético.
Según el comando el mayor indice de parentezco aceptado sería <0.2


## Unir los datos locales con 1000G

Elimine variantes duplicadas del set ChileGenomico

```
plink --bfile chilean_all48_hg19_10 --list-duplicate-vars suppress-first
31640 MB RAM detected; reserving 15820 MB for main workspace.
458097 variants loaded from .bim file.
42 people (24 males, 18 females) loaded from .fam.
42 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 42 founders and 0 nonfounders present.
Calculating allele frequencies... done.
458097 variants and 42 people pass filters and QC.
Among remaining phenotypes, 0 are cases and 42 are controls.
--list-duplicate-vars report written to plink.dupvar .
```

```
plink --bfile chilean_all48_hg19_10 --exclude plink.dupvar --make-bed --out chilean_all48_hg19_11

31640 MB RAM detected; reserving 15820 MB for main workspace.
458097 variants loaded from .bim file.
42 people (24 males, 18 females) loaded from .fam.
42 phenotype values loaded from .fam.
--exclude: 458044 variants remaining.
Warning: At least 1 duplicate ID in --exclude file.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 42 founders and 0 nonfounders present.
Calculating allele frequencies... done.
458044 variants and 42 people pass filters and QC.
Among remaining phenotypes, 0 are cases and 42 are controls.
--make-bed to chilean_all48_hg19_11.bed + chilean_all48_hg19_11.bim +
chilean_all48_hg19_11.fam ... done.
```

Extraiga las variantes presentes en los datos de ChileGenomico.

```
cut -f 2 chilean_all48_hg19_11.bim | sort -u > chilean_all48_hg19_11.snps
```

Extraiga las variantes presentes en los datos de 1000G.

```
cut -f 2 $G/1kG_MDS5.bim | sort -u  > 1kG_MDS5.snps
```

Encuentre la lista de SNPs en común entre ambos sets de datos.

```
comm -12 chilean_all48_hg19_11.snps 1kG_MDS5.snps > common_snps.txt
```

Extraiga los SNPs en común del set 1000G

```
plink --bfile $G/1kG_MDS5 --extract common_snps.txt --recode --make-bed --out 1kG_MDS6
31640 MB RAM detected; reserving 15820 MB for main workspace.
5808310 variants loaded from .bim file.
629 people (0 males, 0 females, 629 ambiguous) loaded from .fam.
Ambiguous sex IDs written to 1kG_MDS6.nosex .
--extract: 314412 variants remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 629 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.999391.
314412 variants and 629 people pass filters and QC.
Note: No phenotypes present.
--make-bed to 1kG_MDS6.bed + 1kG_MDS6.bim + 1kG_MDS6.fam ... done.
--recode ped to 1kG_MDS6.ped + 1kG_MDS6.map ... done.
```

Extraiga los SNPs en común del set ChileGenomico
```
plink --bfile chilean_all48_hg19_11 --extract common_snps.txt --make-bed --out chilean_all48_hg19_12
31640 MB RAM detected; reserving 15820 MB for main workspace.
458044 variants loaded from .bim file.
42 people (24 males, 18 females) loaded from .fam.
42 phenotype values loaded from .fam.
--extract: 314412 variants remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 42 founders and 0 nonfounders present.
Calculating allele frequencies... done.
314412 variants and 42 people pass filters and QC.
Among remaining phenotypes, 0 are cases and 42 are controls.
--make-bed to chilean_all48_hg19_12.bed + chilean_all48_hg19_12.bim +
chilean_all48_hg19_12.fam ... done.
```

## Homologar versión del genoma
```
awk '{print $2, $4}' chilean_all48_hg19_12.bim> buildhapmap.txt
```
buildhapmap.txt contiene un ID y una posición física por SNP en cada línea.

```
plink --bfile 1kG_MDS6 --update-map buildhapmap.txt --make-bed --out 1kG_MDS7
31640 MB RAM detected; reserving 15820 MB for main workspace.
314412 variants loaded from .bim file.
629 people (0 males, 0 females, 629 ambiguous) loaded from .fam.
Ambiguous sex IDs written to 1kG_MDS7.nosex .
--update-map: 314412 values updated.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 629 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.999391.
314412 variants and 629 people pass filters and QC.
Note: No phenotypes present.
--make-bed to 1kG_MDS7.bed + 1kG_MDS7.bim + 1kG_MDS7.fam ... done.
bioinfo1@genoma:~/jtorresnunez/Unidad2/Sesion2$ 
```

## Fusionar los conjuntos de datos HapMap y 1000 Genomes

- Establecer el genoma de referencia en ChileGenomico usando 1000G:
```
awk '{print $2, $5}' 1kG_MDS7.bim> 1kG_ref-list.txt
plink --bfile chilean_all48_hg19_12 --reference-allele 1kG_ref-list.txt --make-bed --out chilean_all48_hg19_13
31640 MB RAM detected; reserving 15820 MB for main workspace.
314412 variants loaded from .bim file.
42 people (24 males, 18 females) loaded from .fam.
42 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 42 founders and 0 nonfounders present.
Calculating allele frequencies... done.
--a1-allele: 314412 assignments made.
314412 variants and 42 people pass filters and QC.
Among remaining phenotypes, 0 are cases and 42 are controls.
--make-bed to chilean_all48_hg19_13.bed + chilean_all48_hg19_13.bim +
chilean_all48_hg19_13.fam ... done.
```

- Resolver problemas de hebra. Compruebe si hay posibles problemas de filamento.

```
awk '{print $2, $5, $6}' 1kG_MDS7.bim> 1kG_MDS7_tmp
awk '{print $2, $5, $6}' chilean_all48_hg19_13.bim > chilean_all48_hg19_13_tmp
sort 1kG_MDS7_tmp chilean_all48_hg19_13_tmp | uniq -u > all_differences.txt
wc -l all_differences.txt
```

## Intercambiar alelos de SNPs con potenciales problemas de hebra

*No aplica para este caso porque no hay problemas de hebra
```
awk '{print $ 1}' all_differences.txt | sort -u > flip_list.txt
plink --bfile chilean_all48_hg19_13 --flip flip_list.txt --reference-allele 1kG_ref-list.txt --make-bed --out chilean_all48_hg19_14
awk '{print $2, $5, $6}' chilean_all48_hg19_14.bim > chilean_all48_hg19_14_tmp
sort 1kG_MDS7_tmp chilean_all48_hg19_14_tmp | uniq -u> uncorresponding_SNPs.txt
wc -l uncorresponding_SNPs.txt
```
```
awk '{print $ 1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exlusion.txt
plink --bfile 1kG_MDS7 --exclude SNPs_for_exlusion.txt --make-bed --out 1kG_MDS8
plink --bfile chilean_all48_hg19_14 --exclude SNPs_for_exlusion.txt --make-bed --out chilean_all48_hg19_15
```

## Combinar 1000G con ChileGenomico.

```
plink --bfile 1kG_MDS8 --bmerge chilean_all48_hg19_15.bed chilean_all48_hg19_15.bim chilean_all48_hg19_15.fam --allow-no-sex --make-bed --out MDS_merge

31640 MB RAM detected; reserving 15820 MB for main workspace.
629 people loaded from 1kG_MDS8.fam.
42 people to be merged from chilean_all48_hg19_15.fam.
Of these, 42 are new, while 0 are present in the base dataset.
314412 markers loaded from 1kG_MDS8.bim.
314412 markers to be merged from chilean_all48_hg19_15.bim.
Of these, 0 are new, while 314412 are present in the base dataset.
Performing single-pass merge (671 people, 314412 variants).
Merged fileset written to MDS_merge-merge.bed + MDS_merge-merge.bim +
MDS_merge-merge.fam .
314412 variants loaded from .bim file.
671 people (24 males, 18 females, 629 ambiguous) loaded from .fam.
Ambiguous sex IDs written to MDS_merge.nosex .
42 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 671 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.999429.
314412 variants and 671 people pass filters and QC.
Among remaining phenotypes, 0 are cases and 42 are controls.  (629 phenotypes
are missing.)
--make-bed to MDS_merge.bed + MDS_merge.bim + MDS_merge.fam ... done.
```

## Analisis de estructura poblacional

- Realizar MDS en datos HapMap-ChileGenomico
```
plink --bfile MDS_merge --extract indepSNP.prune.in --genome --out MDS_merge2

31640 MB RAM detected; reserving 15820 MB for main workspace.
314412 variants loaded from .bim file.
671 people (24 males, 18 females, 629 ambiguous) loaded from .fam.
Ambiguous sex IDs written to MDS_merge2.nosex .
42 phenotype values loaded from .fam.
--extract: 70534 variants remaining.
Using up to 15 threads (change this with --threads).
Before main variant filters, 671 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.99949.
70534 variants and 671 people pass filters and QC.
Among remaining phenotypes, 0 are cases and 42 are controls.  (629 phenotypes
are missing.)
IBD calculations complete.  
Finished writing MDS_merge2.genome .
```

```
plink --bfile MDS_merge --read-genome MDS_merge2.genome --cluster --mds-plot 10 --out MDS_merge2

31640 MB RAM detected; reserving 15820 MB for main workspace.
314412 variants loaded from .bim file.
671 people (24 males, 18 females, 629 ambiguous) loaded from .fam.
Ambiguous sex IDs written to MDS_merge2.nosex .
42 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 671 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.999429.
314412 variants and 671 people pass filters and QC.
Among remaining phenotypes, 0 are cases and 42 are controls.  (629 phenotypes
are missing.)
Clustering... done.                        
Cluster solution written to MDS_merge2.cluster1 , MDS_merge2.cluster2 , and
MDS_merge2.cluster3 .
Performing multidimensional scaling analysis (SVD algorithm, 10
dimensions)... done.
MDS solution written to MDS_merge2.mds .
```

- Generar un archivo con información de poblaciones

```
awk '{print$1,$1,$2}' $G/20100804.ALL.panel > ethnicity_1kG.txt
sed 's/JPT/ASN/g' ethnicity_1kG.txt>ethnicity_1kG2.txt
sed 's/ASW/AFR/g' ethnicity_1kG2.txt>ethnicity_1kG3.txt
sed 's/CEU/EUR/g' ethnicity_1kG3.txt>ethnicity_1kG4.txt
sed 's/CHB/ASN/g' ethnicity_1kG4.txt>ethnicity_1kG5.txt
sed 's/CHD/ASN/g' ethnicity_1kG5.txt>ethnicity_1kG6.txt
sed 's/YRI/AFR/g' ethnicity_1kG6.txt>ethnicity_1kG7.txt
sed 's/LWK/AFR/g' ethnicity_1kG7.txt>ethnicity_1kG8.txt
sed 's/TSI/EUR/g' ethnicity_1kG8.txt>ethnicity_1kG9.txt
sed 's/MXL/AMR/g' ethnicity_1kG9.txt>ethnicity_1kG10.txt
sed 's/GBR/EUR/g' ethnicity_1kG10.txt>ethnicity_1kG11.txt
sed 's/FIN/EUR/g' ethnicity_1kG11.txt>ethnicity_1kG12.txt
sed 's/CHS/ASN/g' ethnicity_1kG12.txt>ethnicity_1kG13.txt
sed 's/PUR/AMR/g' ethnicity_1kG13.txt>ethnicity_1kG14.txt

awk '{if($1~/CDSJ/) pop="MAP"}{if($1~/ARI/) pop="AYM"} {print $1, $2, pop}' chilean_all48_hg19_14.fam > ethnicityfile_CLG.txt

cat ethnicity_1kG14.txt ethnicityfile_CLG.txt | sed -e '1i \ FID IID ethnicity'> ethnicityfile.txt
```

- Graficar resultados de MDS

```
Rscript $W/MDS_merged.R
 IID  FID.x SOL         C1          C2        C3           C4           C5
1 ARI002 ARI002   0 -0.0468148 -0.00396006 0.0793811  0.000967135 -0.001602590
2 ARI003 ARI003   0 -0.0511959 -0.01336250 0.0902059  0.000159820  0.000998344
3 ARI004 ARI004   0 -0.0611720 -0.02932550 0.1097400 -0.003464750 -0.002269980
4 ARI006 ARI006   0 -0.0510816 -0.02300040 0.1025010 -0.002027410 -0.002373030
5 ARI008 ARI008   0 -0.0519357 -0.01140780 0.0927078  0.001383930 -0.000461011
6 ARI010 ARI010   0 -0.0453211  0.01213030 0.0637703  0.004412130 -0.000770667
            C6           C7           C8          C9          C10  FID.y
1  0.002183410 -2.02577e-03  0.000259029 -0.00437632  0.004357640 ARI002
2  0.000106360  6.27563e-05 -0.003669320 -0.00375669 -0.000523978 ARI003
3 -0.000036067 -1.33586e-04  0.001845780  0.00071749 -0.000801809 ARI004
4  0.001486100  3.56327e-03 -0.000320211 -0.00298926 -0.002142250 ARI006
5  0.001165280 -1.87450e-03 -0.000951174 -0.00477461 -0.002649800 ARI008
6 -0.003329840 -1.41725e-02  0.001394420  0.00653913 -0.005424780 ARI010
  ethnicity
1       AYM
2       AYM
3       AYM
4       AYM
5       AYM
6       AYM
```

#### Es la tabla de datos resultante de unir los resultados de un Análisis de Componentes Principales (PCA) (el archivo .eigenvec generado por PLINK) con la información poblacional o metadatos de tus muestras en R.

- Para visualizar PCA creado: 
```
okular MDS.pdf &
```
TAREA
En R, genere gráficos similares para las combinaciones Component 2 vs 3 y 3 vs 4. ¿Qué puede concluir de estos gráficos?

Se exporto archivos MDS_merge2.mds y ethnicityfile.txt, en conjunto con la modificación en el script MDS_merged.r en carpeta ~/WK2019-PopGeno.

[Ver gráfico MDS Componente 2 vs 3](MDS_2vs3.pdf)

Al analizar este grafico destaca una similitud de componentes Aymaras y Mapuches, separados de componentes europeos y asiáticos. 

[Ver gráfico MDS Componente 3 vs 4](MDS_3vs4.pdf)

Al analziar este grafico se evidencia la separaión de 2 subgrupos europeos, lo que podria ser explciado por diferencias geográficas de los individuos.


## Realizar un análisis de ancestría

```
plink --bfile MDS_merge --extract indepSNP.prune.in --make-bed --out MDS_merge_r2_lt_0.2
```

1. ¿Cuántos SNPs quedaron luego del filtro?
quedaron 70534 SNPs

2. ADMIXTURE asume que los individuos no están emparentados. Sin embargo, no realizamos ningún filtro. ¿Por qué?
Debido a que el filtro fue realizado en etapas previas, al evaluar el equilobrio de desligamiento (LD)

```
module load admixture
admixture -j4 --cv MDS_merge_r2_lt_0.2.bed 3 > MDS_merge_r2_lt_0.2.K3.log
for k in $(seq 4 6)
  do
    admixture -j44 --cv MDS_merge_r2_lt_0.2.bed $k > MDS_merge_r2_lt_0.2.K$k.log        
  done
```



Conclusiones

En este trabajo se ejecutó un pipeline bioinformático reproducible de control de calidad y análisis de estructura genómica, aplicando filtros de tasa de datos perdidos, comprobación de sexo, equilibrio de Hardy-Weinberg y remoción de parentesco para depurar la matriz de datos. Esto subraya la relevancia de un estricto filtrado previo para garantizar resultados confiables. Al analizar los gráficos de componentes principales (MDS), se evidenciaron las distintas ancestrías continentales y cómo, al alternar los componentes evaluados, se revelan subgrupos específicos y grados de cercanía genómica entre las poblaciones. Finalmente, la ejecución de ADMIXTURE no pudo concretarse debido a, según lo revisado en web y foros, una posible limitación de memoria RAM en el servidor, un aspecto técnico que se mantendrá en observación para corregir la asignación de recursos en las siguientes sesiones.
