# Tarea sesión 2.2
Bioinformática y análisis genómicos reproducibles

Joaquín Torres Nuñez

Profesor: Ricardo Verdugo

09/10/2026
---
Ejercicios

Consulta el manual de VCFtools y escribe un script que responda lo siguiente para el archivo GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf sin copiarlo a su directorio:

1. ¿Cuántos individuos y variantes (SNPs) tiene el archivo?
18 Individuos
4450360 Variantes (SNPS)
Coamndo utilizado:
```
vcftools --vcf /datos/compartido/ChileGenomico/GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf
```
Analiza el archivo sin modificarlo, lee la cabecera y los registros del archivo VCF desde su ruta original sin realizar cambios. Cuenta el total de individuos (muestras) y el total de variantes (sitios/SNPs) presentes.

2. ¿Cuántos sitios del archivo no tienen datos perdidos?
Se usa comando --max-missing 1 <float>, la cual permite saber aquellas variantes con un call rate de 100%
```
vcftools --vcf /datos/compartido/ChileGenomico/GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf --max-missing 1
```
Despues del filtro,  382626 sitios del archivo no tienen datos perdidos

3. Genera un archivo en tu carpeta de trabajo Sesion2/VCF/data que contenga solo SNPs en una ventana de 2Mb en cualquier cromosoma. Nombra el archivoCLG_Chr<X>_<Start>-<End>Mb.vcf donde es número del cromosoma, es el inicio de la ventana genómica y es el final en megabases.
```
mkdir -p Sesion2/VCF/data
vcftools --vcf /datos/compartido/ChileGenomico/GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf \
  --chr 1 --from-bp 1000000 --to-bp 3000000 --remove-indels \
  --recode --recode-INFO-all --out Sesion2/VCF/data/CLG_Chr1_1-3Mb
mv Sesion2/VCF/data/CLG_Chr1_1-3Mb.recode.vcf Sesion2/VCF/data/CLG_Chr1_1-3Mb.vcf
```
Comandos:
- --chr 1 --from-bp 1000000 --to-bp 3000000 Permite obtener ventana de 2 mb
- --remove-indels Elimina inserciones y deleciones
- --recode --recode-INFO-all Se genera nuevo archivo vcf con los SNPS filtrados
- --out Sesion2/VCF/data/CLG_Chr1_1-3Mb Da la ruta al nuevo archivo
4. Reporta cuántas variantes tienen el archivo generado
```
vcftools --vcf Sesion2/VCF/data/CLG_Chr1_1-3Mb.vcf
```
After filtering, kept 2378 out of a possible 2378 Sites

5. Reporta la cobertura promedio para todos los individuos del set de datos
```
vcftools --vcf Sesion2/VCF/data/CLG_Chr1_1-3Mb.vcf --depth --out Sesion2/VCF/data/cobertura
cat Sesion2/VCF/data/cobertura.idepth
```
Cobertura promedio: 2.03928

6. Calcula la frecuencia de cada alelo para todos los individuos dentro del archivo y guarda el resultado en un archivo
Se utiliza para esta acción el comando --freq. Este comando entrega la frecuencia de alelos para cada sitio en un archivo .frq.
```
vcftools --vcf Sesion2/VCF/data/CLG_Chr1_1-3Mb.vcf --freq --out Sesion2/VCF/data/frecuenciasalelicas
```
Para este caso se crea un archivo con el nombre: frecuenciasalelicas.frq

7. Filtra el archivo de frecuencias para solo incluir variantes bialélicas (tip: awk puede ser útil para realizar esta tarea, tip2: puedes usar bcftools para filtrar variantes con más de dos alelos antes de calcular las frecuencias)
```
vcftools --vcf CLG_Chr1_1-3Mb.vcf --min-alleles 2 --max-alleles 2 --freq --out frecuencias_bialelicas
```
es útil para realizar esta acción los comandos --min-alleles <integer> --max-alleles <integer>. Inclye solo sitios con alelos minimos y maximos que uno defina, para este caso al ser solo bialelicas el minimo y maximo es 2. 
Respuesta: 2377 variantes son bialelicas.

8. Llama a un script escrito en lenguaje R que lee el archivo de frecuencias de variantes bialélicas y guarda un histograma con el espectro de MAF para las variantes bialélicas

9. ¿Cuántos sitios tienen una frecuencia del alelo menor <0.05?
Para realizar esta acción es util el comando --max-maf <float>
```
vcftools --vcf CLG_Chr1_1-3Mb.vcf --min-alleles 2 --max-alleles 2 --max-maf 0.04999 --out maf_menor_005
```
Dicho comando incluye únicamente los sitios con una Frecuencia del Alelo Menor (MAF) mayor o igual al valor de --maf y menor o igual al valor de --max-maf.
Respuesta: 315 sitios tienen un frecuencia alelica menor a 0.05.

10. Calcula la heterocigosidad de cada individuo.
Para realizar dicha acción es util el comando --het el cual calcula la heterocigocidad para cada individuo.
```
vcftools --vcf CLG_Chr1_1-3Mb.vcf --het --out heterocigosidad_individuos
```
11. Calcula la diversidad nucleotídica por sitio.
```
vcftools --vcf CLG_Chr1_1-3Mb.vcf --site-pi --out diversidad_sitio
```

12. 
```
vcftools --vcf CLG_Chr1_1-3Mb.vcf --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --recode-INFO-all --out CLG_Chr1_1-3Mb_maffiltradas
```
Se crea un archivo CLG_Chr1_1-3Mb_maffiltradas con 2062 variantes post-filtracion MAF <0.05.

13. Convierte el archivo GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf a formato plink.
```
plink --vcf GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf \
      --make-bed \
      --out GATK_ChGdb_recalibrated.autosomes.12262013.snps.known
```

Fin parte **vcftools**
---

1. Enlista los archivos plink que hay en `data`.  ¿Qué tipos de archivos son cada uno?
```
mkdir -p data
cp /datos/compartido/ChileGenomico/chilean_all48_hg19* data/
```
chilean_all48_hg19.bed ; archivo .bed, binadio de genotipos
chilean_all48_hg19.bim; archivo .bim archivo de tecto plano con SNPS
chilean_all48_hg19.fam; archivo .fam archivo de texto plano con la información de parentezco entre individuos
chilean_all48_hg19_popinfo.csv; archivo .csv

2. Consulta el manual de [plink1.9](https://www.cog-genomics.org/plink/1.9/formats) y contesta utilizando comandos de plink lo siguiente. Deposita cualquier archivo que generes an una carpeta `Unididad2/Sesion2/PINK/results`:
a. Transforma de formato bed a formato ped (pista: sección Data Managment). El nombre del output debe ser igual, solo cambiando la extensión.
```
plink --bfile chilean_all48_hg19 --allow-extra-chr --recode --out results/chilean_all48_hg19
```


b. Crea otro archivo ped (ojo PPPPed) pero esta vez filtrando los SNPs cuya frecuencia del alelo menor sea menor a 0.05 Y filtrando los individuos con más de 10% missing data. Tu output debe llamarse `chilean_all48_hg19_maf05_missing10`

¿Cuántos SNPs y cuántos individuos fueron removidos por los filtros?
```
plink --bfile chilean_all48_hg19 --allow-extra-chr --maf 0.05 --mind 0.1 --recode --out results/chilean_all48_hg19_maf05_missing10
```
Respuesta: 2 individuos removidos, 347070 variantes removidas.


c. Realiza un reporte de equilibrio de Hardy-Weinberg sobre el archivo `chilean_all48_hg19_maf05_missing10` creado en el ejercicio anterior. El nombre del archivo de tu output debe contener chilean_all48_hg19_maf05_missing10.
```
plink --file results/chilean_all48_hg19_maf05_missing10 --allow-extra-chr --hardy --out results/chilean_all48_hg19_maf05_missing10
```
Cda columna corresponde a:
-CHR: Cromosoma
- SNP: identificador SNPS
- Test: indica que subgrupo fue analizado
- A1/A2: Alelos
- Geno: Numero de genotipos
- O(HET): eheterocigosidad observada
- E(HET): Heterocigosidad esperada
- P: Valor de Hardy Weinberg

d. Observa el archivo `chilean_all48_hg19.fam`. Consulta la documentación de plink para determinar que es cada columna. ¿Qué información hay y no hay en este archivo?

```
head chilean_all48_hg19.fam
```
Las columnas corresponden a:
-FID: Family ID
- IID: Individual ID
-PID: ID padre
-MID: ID madre
-Sexo; Columna 5
- Fenotipo: columna 6

4. Utiliza la info el archivo `data/chilean_all48_hg19_popinfo.csv` y el comando `update-ids` de plink para cambiar los nombres de las muestras de `data/chilean_all48_hg19.fam` de tal forma que el family ID corresponda a la info de la columna `Ancestry` en `chilean_all48_hg19_popinfo.csv`. Pista: este ejercicio requiere varias operaciones, puedes dividirlas en diferentes scripts de bash o de R y bash. Tu respuesta debe incluir todos los scripts (y deben estar en /code). 

5) Realiza un cuna comparación entre el sexo y archivo `fam`y el `popinfo` y calcula la proporción de discordancias
6) Realiza un test de estimación de sexo usando plink y reporta los resultados en formato de tabla para todos los individuos con discordancia entre el sexto reportado en `fam` y el calculado con plink.
7) Genera una tabla de contingencia de individuos por sexo y ancestría (hint: ver columna Ancestry en el archivo `popinfo`)




## Conclusión

En el presente trabajo se realizo ejercicios y analisis de secuencias utilizando VCFtools y PLINK. VCFtools corresponde a una herramienta de línea de comandos cuyo fin es trabajar archivos con formato VCF, los cuales son utilizados en forma amplia para almacenar SNPS genómicos. A través de VCFtools se evaluaron métricas clave como la cobertura promedio, la diversidad nucleotídica y la distribución de frecuencias alélicas (MAF), aislando variantes bialélicas y sitios sin datos perdidos. Destaco la facilidad de la ejecución de comandos y buena explicación del manual VCFtools que facilito esta parte de la tarea. Con PLINK, herramienta utilizada para el análisis de genética de poblaciones, se gestionó la transición entre archivos binarios (.bed/.bim/.fam) y de texto (.ped/.map), aplicando filtros rigurosos (MAF ≥ 0.05 y periddas ≤ 10%) que removieron 2 individuos y 347,070 variantes raras o de baja calidad previo a la evaluación del equilibrio de Hardy-Weinberg. La principal dificultad del práctico radicó en la comprensión de la estructura de los archivos generados, lo que conllevó un mayor tiempo de análisis e impidió completar las últimas actividades del tutorial. Como aprendizaje principal, se evidencia la necesidad de dominar exhaustivamente la lógica del pipeline y los formatos de entrada y salida en genética de poblaciones antes de ejecutar los comandos, asegurando así una correcta interpretación de los resultados.



