# Population Genetics Practical


Bioinformática y análisis genómicos reproducibles
Joaquín Torres N.
Profesores:
- Constanza de la Fuente
- Ricardo Verdugo

09/11/2026

---
## PROJECT 1: Do all present-day populations from Europe display the same 3-way admixture?

Pregunta de investigación: ¿Todas las poblaciones europeas modernas derivan del mismo mestizaje de 3 vías?

Objetivo: Evaluar si las poblaciones europeas presentes en el estudio poseen la misma mixtura de 3 vías. Se analizaron las poblaciones de Gran Bretaña, Italia, Finlandia, Iberica-mediterrena y Sardinia.

Metodos: La metodología empleada en este proyecto corresponde al análisis estadístico de genómica poblacional 
y ancestría mediante el uso del paquete admixtools en R, fundamentado en los estadísticos de frecuencias alélicas (estadísticos f2, f3 y f4)

1. Cargar librerias
```
library(admixtools)
install.packages("dplyr")
library(dplyr)
```

Debido a problemas con la instalación de tidyverse se decide solo instalar paquete **dplyr**

2. Cargar f2 blocks
Durante la sesión presencial fueron calculados f2. Se agregaron a directorio de trabajo. 

```
library(admixtools)
outdir <- "/home/joaquin/Escritorio/genantiguo/aadr_1000G_f2_all"
f2_blocks <- f2_from_precomp(outdir)
```

3. Definir poblaciones antiguas de referencia. A partir de estas poblaciones se modelará la población target.  
```
source1 <- c("Turkey_Marmara_Barcin_N.AG", 
             "Russia_Samara_EBA_Yamnaya.AG", 
             "Luxembourg_Mesolithic.DG")

```

4. Definir outgroups. Se designan 5 grupos outgroups para el calculo de deriva genética. Permite la cuantificación del flujo de ancestría. Se calcula la deriva genetica entre grupos desde la separación con el outgroup.

```
outgroup1 <- c("Mbuti.DG", "CHB.DG", "Papuan.DG", "Russia_UstIshim_IUP.DG", "Denisova.DG")
```

5. Concatenación de poblaciones europeas para comparar. 

targets <- c("GBR.DG", "IBS.DG", "TSI.DG", "FIN.DG", "Sardinian.DG")

6. Calculo de f3. Se utiliza admixtools para el calculo de la deriva genetica compartida entre la población de estudio (GBR.DG, IBS.DG, TSI.DG, FIN.DG, Sardinian.DG) con la población de comparación (Turkey_Marmara_Barcin_N.AG Russia_Samara_EBA_Yamnaya.AG, Luxembourg_Mesolithic.DG) utilizando como outgroup la población Mbuti. 
```
f3_results <- f3(f2_blocks, pop1 = "Mbuti.DG", pop2 = "GBR.DG", pop3 = source1)
```
**Resultado cálculo f3:**
```
> print(f3_results)
A tibble: 3 × 7
  pop1     pop2   pop3                            est       se     z     p
  <chr>    <chr>  <chr>                         <dbl>    <dbl> <dbl> <dbl>
1 Mbuti.DG GBR.DG Luxembourg_Mesolithic.DG     0.0745 0.000553  135.     0
2 Mbuti.DG GBR.DG Russia_Samara_EBA_Yamnaya.AG 0.0719 0.000493  146.     0
3 Mbuti.DG GBR.DG Turkey_Marmara_Barcin_N.AG   0.0714 0.000487  147.     0
```
Del presente resultado destacan los estadigrafos est (calculo f3) y Z score. Se observa que la deriva génica es mayor entre las población GRB y luxemburgo mesolitica, en comparacion a las poblaciones Yamnaya y Barcin, no obstante sugiere una mezcla de las 3 poblaciones en el genoma de GBR. Al analizar el estadígrafo Z-score se aprecia que los valores observados poseen una alta confiabilidad, confirmando estadisticamente el flujo genético de las  poblaciones source a la población target. 

7. Calculo de f4 o ABBA BABA. Detecta el flujo génico entre 4 poblaciones. Se agrega en una rama la presencia de un outgroup (Mbuti.DG, pop4)
```
f4_results <- f4(f2_blocks, 
                 pop1 = "GBR.DG", 
                 pop2 = ref_pops, 
                 pop3 = ref_pops, 
                 pop4 = "Mbuti.DG")
```

**Resultado:**
```
  pop1   pop2                         pop3                         pop4           est       se      z         p
  <chr>  <chr>                        <chr>                        <chr>        <dbl>    <dbl>  <dbl>     <dbl>
1 GBR.DG Luxembourg_Mesolithic.DG     Luxembourg_Mesolithic.DG     Mbuti.DG -0.0311   0.00112  -27.9  6.49e-171
2 GBR.DG Luxembourg_Mesolithic.DG     Turkey_Marmara_Barcin_N.AG   Mbuti.DG  0.00156  0.000305   5.13 2.91e-  7
3 GBR.DG Luxembourg_Mesolithic.DG     Russia_Samara_EBA_Yamnaya.AG Mbuti.DG -0.000575 0.000324  -1.77 7.63e-  2
4 GBR.DG Turkey_Marmara_Barcin_N.AG   Luxembourg_Mesolithic.DG     Mbuti.DG  0.00461  0.000198  23.4  1.29e-120
5 GBR.DG Turkey_Marmara_Barcin_N.AG   Turkey_Marmara_Barcin_N.AG   Mbuti.DG -0.00684  0.000145 -47.1  0        
6 GBR.DG Turkey_Marmara_Barcin_N.AG   Russia_Samara_EBA_Yamnaya.AG Mbuti.DG  0.00317  0.000135  23.5  9.33e-122
7 GBR.DG Russia_Samara_EBA_Yamnaya.AG Luxembourg_Mesolithic.DG     Mbuti.DG  0.00198  0.000241   8.20 2.42e- 16
8 GBR.DG Russia_Samara_EBA_Yamnaya.AG Turkey_Marmara_Barcin_N.AG   Mbuti.DG  0.00267  0.000150  17.8  8.40e- 71
9 GBR.DG Russia_Samara_EBA_Yamnaya.AG Russia_Samara_EBA_Yamnaya.AG Mbuti.DG -0.00657  0.000196 -33.5  5.58e-246
```
Debido a los resultados observados se observa en la fila 2 un valor z=5.13, esto orienta una mayor cantidad de alelos comaprtidos entre la población GBR y Barcín. No ibstante al evaluar la fila 4 se observa tambien resultado Z significativo entre GBR y población luxemburgo mesolitica (Z-score 23.4). Esto sustenta el resultado anterior en donde se observa y sugiere una poblacion GBR mezclada por estos dos grupos. Finalmente al analizar la fila 6 se observa que la poblacion GBR tambien tiene una fuerte relación con la poblacipon Yamnaya, pero con un menor numero de alelos compartidos en comapración a las poblaciones Mesoliticas de Luxemburgo y población Barcín.


8. Este comando ejecuta el análisis qpWave, una prueba estadística para determinar el número mínimo de corrientes de ancestría independientes necesarias para explicar la relación entre un grupo de poblaciones de interés (left) y un outgroups (right). 

```
wave1 <- qpwave(f2_blocks,
                left = c("GBR.DG", source1),
                right = outgroup1)
```
**Resultado:**

```
$rankdrop
# A tibble: 3 × 7
  f4rank   dof   chisq        p dofdiff chisqdiff  p_nested
   <int> <int>   <dbl>    <dbl>   <int>     <dbl>     <dbl>
1      2     2   0.354 8.38e- 1       4      38.6  8.38e- 8
2      1     6  39.0   7.27e- 7       6     292.   4.43e-60
3      0    12 331.    1.51e-63      NA      NA   NA 
```      
Se evalua el numero mínimo de ancestrias necesarias para explciar la relacion genética entre el grupo target (left) con las poblaciones a comparar (right). Se observa que la presencia de 3 (p=0.838) flujos ancestrales son suficientes para explicar la similitud alelica (k+1). Esto confirma el mestizaje de 3 vías presente en la población de Gran Bretaña.

9. Este código ejecuta dos modelos de mezcla genética (admixture modeling) con qpAdm para estimar de forma cuantitativa las proporciones de ancestría del grupo británico moderno (GBR.DG), comparando una hipótesis de 2 fuentes contra una de 3 fuentes.

```
admix_2way <- qpadm(f2_blocks, 
                    left = source1[1:2], 
                    right = outgroup1, 
                    target = "GBR.DG")

admix_3way <- qpadm(f2_blocks, 
                    left = source1, 
                    right = outgroup1, 
                    target = "GBR.DG")

```
**Resultado Modelo 2-way:**
```
> print(admix_2way$weights)
# A tibble: 2 × 5
  target left                         weight     se     z
  <chr>  <chr>                         <dbl>  <dbl> <dbl>
1 GBR.DG Turkey_Marmara_Barcin_N.AG    0.280 0.0496  5.63
2 GBR.DG Russia_Samara_EBA_Yamnaya.AG  0.720 0.0496 14.5 
```

**Resultado Modelo 3-way:**
```
> print(admix_3way$weights)
# A tibble: 3 × 5
  target left                         weight     se     z
  <chr>  <chr>                         <dbl>  <dbl> <dbl>
1 GBR.DG Turkey_Marmara_Barcin_N.AG    0.290 0.0572  5.08
2 GBR.DG Russia_Samara_EBA_Yamnaya.AG  0.846 0.108   7.81
3 GBR.DG Luxembourg_Mesolithic.DG     -0.137 0.0921 -1.48
```

AL analizar el estudio de mestizaje destaca un peso negativo para la población Mesolítica de Luxemburgo. Dicho evento se puede deber a que los allelos comaprtidos llegaron meidante otras poblaciones a la poblacion GBR. 

10. Posterior a esto se modificó el script para que realice todos los calculos de las poblaciones target a partir de una iteración. Los resultados finales de admixture se muestran a continuación:

10.1 GBR.DG 
```
Modelo 2-way 
# A tibble: 2 × 5
  target left                         weight     se     z
  <chr>  <chr>                         <dbl>  <dbl> <dbl>
1 GBR.DG Turkey_Marmara_Barcin_N.AG    0.280 0.0496  5.63
2 GBR.DG Russia_Samara_EBA_Yamnaya.AG  0.720 0.0496 14.5 
```
```
Modelo 3-way 
# A tibble: 3 × 5
  target left                         weight     se     z
  <chr>  <chr>                         <dbl>  <dbl> <dbl>
1 GBR.DG Turkey_Marmara_Barcin_N.AG    0.290 0.0572  5.08
2 GBR.DG Russia_Samara_EBA_Yamnaya.AG  0.846 0.108   7.81
3 GBR.DG Luxembourg_Mesolithic.DG     -0.137 0.0921 -1.48
```
10.2 IBS.DG
```
Modelo 2-way
  target left                         weight     se     z
  <chr>  <chr>                         <dbl>  <dbl> <dbl>
1 IBS.DG Turkey_Marmara_Barcin_N.AG    0.669 0.0450 14.9 
2 IBS.DG Russia_Samara_EBA_Yamnaya.AG  0.331 0.0450  7.35
```
```
Modelo 3-way 
  target left                         weight     se     z
  <chr>  <chr>                         <dbl>  <dbl> <dbl>
1 IBS.DG Turkey_Marmara_Barcin_N.AG    0.681 0.0609 11.2 
2 IBS.DG Russia_Samara_EBA_Yamnaya.AG  0.618 0.119   5.20
3 IBS.DG Luxembourg_Mesolithic.DG     -0.299 0.106  -2.82
```

10.3 TSI.DG
```
Modelo 2-way 
  target left                         weight     se     z
  <chr>  <chr>                         <dbl>  <dbl> <dbl>
1 TSI.DG Turkey_Marmara_Barcin_N.AG    0.577 0.0502 11.5 
2 TSI.DG Russia_Samara_EBA_Yamnaya.AG  0.423 0.0502  8.43
```
```
Modelo 3-way 
  target left                         weight     se     z
  <chr>  <chr>                         <dbl>  <dbl> <dbl>
1 TSI.DG Turkey_Marmara_Barcin_N.AG    0.609 0.0608 10.0 
2 TSI.DG Russia_Samara_EBA_Yamnaya.AG  0.698 0.121   5.77
3 TSI.DG Luxembourg_Mesolithic.DG     -0.307 0.106  -2.91
```
10.4 FIN.DG
```
Modelo 2-way 
  target left                         weight     se     z
  <chr>  <chr>                         <dbl>  <dbl> <dbl>
1 FIN.DG Turkey_Marmara_Barcin_N.AG   -0.412 0.0911 -4.52
2 FIN.DG Russia_Samara_EBA_Yamnaya.AG  1.41  0.0911 15.5 
```
```
Modelo 3-way
  target left                         weight    se     z
  <chr>  <chr>                         <dbl> <dbl> <dbl>
1 FIN.DG Turkey_Marmara_Barcin_N.AG   -0.344 0.124 -2.77
2 FIN.DG Russia_Samara_EBA_Yamnaya.AG  1.77  0.265  6.68
3 FIN.DG Luxembourg_Mesolithic.DG     -0.430 0.247 -1.74
```

10.5 Sardinian.DG
```
Modelo 2-way 
  target       left                         weight     se     z
  <chr>        <chr>                         <dbl>  <dbl> <dbl>
1 Sardinian.DG Turkey_Marmara_Barcin_N.AG    0.753 0.0437 17.3 
2 Sardinian.DG Russia_Samara_EBA_Yamnaya.AG  0.247 0.0437  5.65
```
```
Modelo 3-way
  target       left                         weight     se     z
  <chr>        <chr>                         <dbl>  <dbl> <dbl>
1 Sardinian.DG Turkey_Marmara_Barcin_N.AG    0.763 0.0580 13.2 
2 Sardinian.DG Russia_Samara_EBA_Yamnaya.AG  0.493 0.116   4.25
3 Sardinian.DG Luxembourg_Mesolithic.DG     -0.256 0.103  -2.50
```


Respuesta: A partir de los resultaods expuestos no toda población europea comparte el mismo linaje ancestral de 3 fuentes. En las poblaciones analizadas se observa que contribuciones diferentes en la mezcla de las poblaciones de referencia, particularmente de la población Luxemburgo mesolitica la cual presento Z scores negativos en los distintos linajes planteados.


## Proyecto 2: Steppe formation model.

Pregunta de investigación: ¿Pueden las poblaciones Steppe Emba (Russia_MLBA_Sintashta, Kazakhstan_Maitan_MLBA_Alakul, Russia_LBA_Srubnaya_Alakul) explicar su mestizaje por un modelamiento de 2 vías entre población Iraní neolítica y Sidelkino? 

Objetivo: Determinar si la población Steppe EMBA puede ser modelada por un mestizaje entre Irani neolitco y población Sidelkino.

Metodología: La metodología empleada en este proyecto corresponde al análisis estadístico de genómica poblacional 
y ancestría mediante el uso del paquete admixtools en R, fundamentado en los estadísticos de frecuencias alélicas (estadísticos f2, f3 y f4)

Fueron calculados de forma previa los f2 de los grupos GanjDareh, Sidelkino (source) y Mitan MLBA (target) dado que no fueron cargados calculados los bloques f2 de estas poblaciones


**Calculo F3**
```

 pop1     pop2                             pop3                      est       se     z     p
  <chr>    <chr>                            <chr>                   <dbl>    <dbl> <dbl> <dbl>
1 Mbuti.DG Kazakhstan_Maitan_MLBA_Alakul.AG Iran_GanjDareh_N.AG    0.0691 0.000645  107.     0
2 Mbuti.DG Kazakhstan_Maitan_MLBA_Alakul.AG Russia_Sidelkino_HG.SG 0.0767 0.000763  101.     0
3 Mbuti.DG Russia_LBA_Srubnaya_Alakul.SG    Iran_GanjDareh_N.AG    0.0687 0.000649  106.     0
4 Mbuti.DG Russia_LBA_Srubnaya_Alakul.SG    Russia_Sidelkino_HG.SG 0.0779 0.000747  104.     0
5 Mbuti.DG Russia_MLBA_Sintashta.AG         Iran_GanjDareh_N.AG    0.0692 0.000627  110.     0
6 Mbuti.DG Russia_MLBA_Sintashta.AG         Russia_Sidelkino_HG.SG 0.0771 0.000736  105.     0
```

Se observa una mayor deriva genica de las poblaciones de la estepa con la población Sidelkino (cazadores-recolectores del este):
En las tres poblaciones de la estepa (Sintashta, Srubnaya, Maitan), la relación con la poblacipon Sidelkino presenta estadisticos f3 mayores (>0.07) en comapracion a la poblacipon Iraní (<0.07). Ambos datos respaldados por un Z-score alto (>100). No obstante se presenta una alta porporción de alelos compartidos de ambas poblacioens de referencia con el grupo target.



**Calculo f4 ABBA BABA**
```

   pop1                             pop2                   pop3                   pop4          est       se       z        p
   <chr>                            <chr>                  <chr>                  <chr>       <dbl>    <dbl>   <dbl>    <dbl>
 1 Russia_MLBA_Sintashta.AG         Iran_GanjDareh_N.AG    Iran_GanjDareh_N.AG    Mbuti.DG -0.0206  0.000504  -40.8  0       
 2 Russia_MLBA_Sintashta.AG         Iran_GanjDareh_N.AG    Russia_Sidelkino_HG.SG Mbuti.DG  0.0101  0.000497   20.4  1.26e-92
 3 Russia_MLBA_Sintashta.AG         Russia_Sidelkino_HG.SG Iran_GanjDareh_N.AG    Mbuti.DG  0.00223 0.000538    4.15 3.31e- 5
 4 Russia_MLBA_Sintashta.AG         Russia_Sidelkino_HG.SG Russia_Sidelkino_HG.SG Mbuti.DG -0.140   0.000610 -230.   0       
 5 Kazakhstan_Maitan_MLBA_Alakul.AG Iran_GanjDareh_N.AG    Iran_GanjDareh_N.AG    Mbuti.DG -0.0207  0.000510  -40.5  0       
 6 Kazakhstan_Maitan_MLBA_Alakul.AG Iran_GanjDareh_N.AG    Russia_Sidelkino_HG.SG Mbuti.DG  0.00973 0.000554   17.6  5.75e-69
 7 Kazakhstan_Maitan_MLBA_Alakul.AG Russia_Sidelkino_HG.SG Iran_GanjDareh_N.AG    Mbuti.DG  0.00214 0.000565    3.79 1.51e- 4
 8 Kazakhstan_Maitan_MLBA_Alakul.AG Russia_Sidelkino_HG.SG Russia_Sidelkino_HG.SG Mbuti.DG -0.141   0.000662 -212.   0       
 9 Russia_LBA_Srubnaya_Alakul.SG    Iran_GanjDareh_N.AG    Iran_GanjDareh_N.AG    Mbuti.DG -0.0211  0.000537  -39.3  0       
10 Russia_LBA_Srubnaya_Alakul.SG    Iran_GanjDareh_N.AG    Russia_Sidelkino_HG.SG Mbuti.DG  0.0109  0.000519   21.1  1.89e-98
11 Russia_LBA_Srubnaya_Alakul.SG    Russia_Sidelkino_HG.SG Iran_GanjDareh_N.AG    Mbuti.DG  0.00172 0.000538    3.19 1.43e- 3
12 Russia_LBA_Srubnaya_Alakul.SG    Russia_Sidelkino_HG.SG Russia_Sidelkino_HG.SG Mbuti.DG -0.139   0.000635 -219.   0       
```

Al evaluar el estadistico f4 se observa valores variados al comparar las poblaciones target con las pobalciones GanjDAreh y Sidelkino. Al evaluar la fila 2 6 y 10 se observa un Z score significativo lo que sugiere una alta deriva genica entre las poblaciones EMPA con la población Sidelkino. Al analizar las filas 3, 7 y 11 se observan de igual forma valores signiicatis entre las pobkaciones target y GanjDareh los Zsocre son de menor valor que los evidenciados anteriormetne, lo que sugiere una menor similitud alelica entre estos dos componentes.

**Calculo qpwave**
```
  f4rank   dof chisq        p dofdiff chisqdiff  p_nested
   <int> <int> <dbl>    <dbl>   <int>     <dbl>     <dbl>
1      1     3  1.77 6.22e- 1       5      84.6  9.01e-17
2      0     8 86.4  2.50e-15      NA      NA   NA       
```

El análisis demuestra de manera formal que las poblaciones estudiadas derivan de 2 flujos genéticos/fuentes ancestrales distintas en relación con el conjunto de outgroups utilizado. Un modelo con una única fuente queda estadísticamente descartado, mientras que un modelo de 2 vías (2-way admixture / 2 streams of ancestry) es estadísticamente suficiente y óptimo.

**Calculo qpadm**
```
  target                           left                   weight     se     z pvalue
  <chr>                            <chr>                   <dbl>  <dbl> <dbl>  <dbl>
1 Russia_MLBA_Sintashta.AG         Iran_GanjDareh_N.AG     0.334 0.0809  4.13  0.622
2 Russia_MLBA_Sintashta.AG         Russia_Sidelkino_HG.SG  0.666 0.0809  8.23  0.622
3 Kazakhstan_Maitan_MLBA_Alakul.AG Iran_GanjDareh_N.AG     0.457 0.0804  5.68  0.900
4 Kazakhstan_Maitan_MLBA_Alakul.AG Russia_Sidelkino_HG.SG  0.543 0.0804  6.75  0.900
5 Russia_LBA_Srubnaya_Alakul.SG    Iran_GanjDareh_N.AG     0.275 0.0860  3.20  0.987
6 Russia_LBA_Srubnaya_Alakul.SG    Russia_Sidelkino_HG.SG  0.725 0.0860  8.43  0.987
```

Las tres poblaciones de la estepa del Bronce Medio/Tardío pueden explicarse con éxito como mezclas genéticas derivadas de ancestros tipo cazadores-recolectores del este (EHG) y de agricultores/poblaciones del suroeste asiático (Neolítico de Irán / CHG-related), pero con valores diferentes, siendo mas predominante de la poblacon Sidelkino que Iraní GanjDareh. 



## Proyecto 3: Peopling of the Americas.

Pregunta de investigación: ¿Como contribuyeron poblaciones americanas antiguas a la estructura genética  poblacioens americanas modernas, y si alguna podría teenr influencia de poblaciones Australo-Melanesia? 

Objetivo: Evaluar modelos de múltiples olas de migración hacia las Américas (por ejemplo, emparentados con Anzick vs. Beringios Antiguos), y determinar si algunos grupos presentan una afinidad adicional con los Australo-Melanesios (una señal controversial).

Metodología: Se calcularon f2, f3, f4, qpwave, qpadm de poblaciones target, referencia y outgroups para este proyecto. Se estudiaron las poblaciones Karitiana, PEL, Pima, Surui y MXL. 

Calculo f3:

**Resultados:**
```
   pop1     pop2         pop3                                est       se     z     p
   <chr>    <chr>        <chr>                             <dbl>    <dbl> <dbl> <dbl>
 1 Mbuti.DG Karitiana.DG USA_Ancient_Beringian.SG         0.0861 0.000703  123.     0
 2 Mbuti.DG Karitiana.DG USA_Anzick_realigned.SG          0.0971 0.000766  127.     0
 3 Mbuti.DG Karitiana.DG USA_Nevada_SpiritCave_11000BP.SG 0.0991 0.000749  132.     0
 4 Mbuti.DG MXL.DG       USA_Ancient_Beringian.SG         0.0735 0.000569  129.     0
 5 Mbuti.DG MXL.DG       USA_Anzick_realigned.SG          0.0788 0.000599  132.     0
 6 Mbuti.DG MXL.DG       USA_Nevada_SpiritCave_11000BP.SG 0.0800 0.000590  136.     0
 7 Mbuti.DG PEL.DG       USA_Ancient_Beringian.SG         0.0810 0.000622  130.     0
 8 Mbuti.DG PEL.DG       USA_Anzick_realigned.SG          0.0899 0.000669  134.     0
 9 Mbuti.DG PEL.DG       USA_Nevada_SpiritCave_11000BP.SG 0.0918 0.000651  141.     0
10 Mbuti.DG Pima.DG      USA_Ancient_Beringian.SG         0.0856 0.000689  124.     0
11 Mbuti.DG Pima.DG      USA_Anzick_realigned.SG          0.0948 0.000746  127.     0
12 Mbuti.DG Pima.DG      USA_Nevada_SpiritCave_11000BP.SG 0.0966 0.000730  132.     0
13 Mbuti.DG Surui.DG     USA_Ancient_Beringian.SG         0.0859 0.000731  117.     0
14 Mbuti.DG Surui.DG     USA_Anzick_realigned.SG          0.0972 0.000797  122.     0
15 Mbuti.DG Surui.DG     USA_Nevada_SpiritCave_11000BP.SG 0.0991 0.000770  129.     0
```

La poblacion Karitiana y Pima comparte mas alelos con la pobalcion ancestral Anzick. La poblacion MXL, Surui y PEL, originarios de México, Brasil y Perú respectivamente poseen mas alelos comaprtidos con la población ancestral Nevada. Esto indica una variación entre las poblacioens indigenas americanas en relación a la estructura de su genoma y posibles ancestrías de estos. 


Calculo f4:
```
**Resultados:**
Relacion con Papuan
 pop1           pop2                             pop3      pop4           est       se      z        p
   <chr>        <chr>                            <chr>     <chr>        <dbl>    <dbl>  <dbl>    <dbl>
 1 Pima.DG      USA_Anzick_realigned.SG          Papuan.DG Mbuti.DG -0.00134  0.000366  -3.66 2.49e- 4
 2 Pima.DG      USA_Ancient_Beringian.SG         Papuan.DG Mbuti.DG -0.00164  0.000395  -4.14 3.40e- 5
 3 Pima.DG      USA_Nevada_SpiritCave_11000BP.SG Papuan.DG Mbuti.DG -0.00154  0.000373  -4.14 3.52e- 5
 4 Karitiana.DG USA_Anzick_realigned.SG          Papuan.DG Mbuti.DG -0.000825 0.000381  -2.17 3.01e- 2
 5 Karitiana.DG USA_Ancient_Beringian.SG         Papuan.DG Mbuti.DG -0.00112  0.000427  -2.63 8.55e- 3
 6 Karitiana.DG USA_Nevada_SpiritCave_11000BP.SG Papuan.DG Mbuti.DG -0.00103  0.000394  -2.60 9.23e- 3
 7 Surui.DG     USA_Anzick_realigned.SG          Papuan.DG Mbuti.DG -0.000546 0.000394  -1.38 1.66e- 1
 8 Surui.DG     USA_Ancient_Beringian.SG         Papuan.DG Mbuti.DG -0.000842 0.000447  -1.88 5.96e- 2
 9 Surui.DG     USA_Nevada_SpiritCave_11000BP.SG Papuan.DG Mbuti.DG -0.000747 0.000385  -1.94 5.23e- 2
10 PEL.DG       USA_Anzick_realigned.SG          Papuan.DG Mbuti.DG -0.00263  0.000338  -7.80 6.30e-15
11 PEL.DG       USA_Ancient_Beringian.SG         Papuan.DG Mbuti.DG -0.00293  0.000373  -7.86 3.93e-15
12 PEL.DG       USA_Nevada_SpiritCave_11000BP.SG Papuan.DG Mbuti.DG -0.00283  0.000333  -8.52 1.59e-17
13 MXL.DG       USA_Anzick_realigned.SG          Papuan.DG Mbuti.DG -0.00565  0.000362 -15.6  5.71e-55
14 MXL.DG       USA_Ancient_Beringian.SG         Papuan.DG Mbuti.DG -0.00595  0.000379 -15.7  1.35e-55
15 MXL.DG       USA_Nevada_SpiritCave_11000BP.SG Papuan.DG Mbuti.DG -0.00585  0.000372 -15.7  9.46e-56
```
**Relacion con fuentes source**
```
pop1            pop2                             pop3                             pop4           est       se        z        p
   <chr>        <chr>                            <chr>                            <chr>        <dbl>    <dbl>    <dbl>    <dbl>
 1 Pima.DG      USA_Anzick_realigned.SG          USA_Anzick_realigned.SG          Mbuti.DG -0.129    0.000558 -231.    0       
 2 Pima.DG      USA_Anzick_realigned.SG          USA_Ancient_Beringian.SG         Mbuti.DG -0.000277 0.000630   -0.440 6.60e- 1
 3 Pima.DG      USA_Anzick_realigned.SG          USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG -0.00246  0.000735   -3.35  8.09e- 4
 4 Pima.DG      USA_Ancient_Beringian.SG         USA_Anzick_realigned.SG          Mbuti.DG  0.00886  0.000735   12.1   1.81e-33
 5 Pima.DG      USA_Ancient_Beringian.SG         USA_Ancient_Beringian.SG         Mbuti.DG -0.122    0.000495 -246.    0       
 6 Pima.DG      USA_Ancient_Beringian.SG         USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG  0.0104   0.000703   14.8   9.64e-50
 7 Pima.DG      USA_Nevada_SpiritCave_11000BP.SG USA_Anzick_realigned.SG          Mbuti.DG -0.00426  0.000705   -6.04  1.50e- 9
 8 Pima.DG      USA_Nevada_SpiritCave_11000BP.SG USA_Ancient_Beringian.SG         Mbuti.DG -0.000513 0.000563   -0.911 3.62e- 1
 9 Pima.DG      USA_Nevada_SpiritCave_11000BP.SG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG -0.112    0.000501 -223.    0       
10 Karitiana.DG USA_Anzick_realigned.SG          USA_Anzick_realigned.SG          Mbuti.DG -0.127    0.000585 -217.    0       
```

Al evaluar en primera instancia una deriva genica proveniente de la población Papuan (Oceanía), esta presenta valores negatiovos al comaprarla con todas las poblaciones target del estudio, lo que indica una ausencia de alelos comaprtidos entre las poblaciones americanas modernas y población oceánica, sugiriendo que las poblaciones americanas ancestrales son la sque poseen un mayor impacto en l estructura genómica de la población moderna. 

Esto se evidencia en la siguiente tabal, donde para ejemplificar se tomará a la pobalción Pima como target. En esta se aprecia una mayor simetría con las poblaciones Anzick y Nevada_Spirit cave en comapración a la poblacipon Beringian.

qpwave

**Resultado:**
```
========================================
Resultados qpWave para: Pima.DG 
========================================
$f4_var_rcond
[1] 0.118186

$f4_var_singular_loadings
NULL

$f4
# A tibble: 9 × 8
  pop1    pop2                             pop3     pop4            est       se     z         p
  <chr>   <chr>                            <chr>    <chr>         <dbl>    <dbl> <dbl>     <dbl>
1 Pima.DG USA_Anzick_realigned.SG          Mbuti.DG CHB.DG      0.00341 0.000344  9.93 3.08e- 23
2 Pima.DG USA_Anzick_realigned.SG          Mbuti.DG Papuan.DG   0.00134 0.000366  3.66 2.49e-  4
3 Pima.DG USA_Anzick_realigned.SG          Mbuti.DG Denisova.DG 0.00948 0.000295 32.1  9.42e-227
4 Pima.DG USA_Ancient_Beringian.SG         Mbuti.DG CHB.DG      0.00332 0.000372  8.91 4.91e- 19
5 Pima.DG USA_Ancient_Beringian.SG         Mbuti.DG Papuan.DG   0.00164 0.000395  4.14 3.40e-  5
6 Pima.DG USA_Ancient_Beringian.SG         Mbuti.DG Denisova.DG 0.00947 0.000323 29.3  1.16e-188
7 Pima.DG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG CHB.DG      0.00434 0.000329 13.2  1.13e- 39
8 Pima.DG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG Papuan.DG   0.00154 0.000373  4.14 3.52e-  5
9 Pima.DG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG Denisova.DG 0.00889 0.000286 31.1  1.38e-212

$rankdrop
# A tibble: 3 × 7
  f4rank   dof    chisq      p dofdiff chisqdiff p_nested
   <int> <int>    <dbl>  <dbl>   <int>     <dbl>    <dbl>
1      2     1    0.490 0.484        3      12.3  0.00637
2      1     4   12.8   0.0123       5    1965.   0      
3      0     9 1978.    0           NA      NA   NA      

```
```
========================================
Resultados qpWave para: Karitiana.DG 
========================================
$f4_var_rcond
[1] 0.09974604

$f4_var_singular_loadings
NULL

$f4
# A tibble: 9 × 8
  pop1         pop2                             pop3     pop4             est       se     z         p
  <chr>        <chr>                            <chr>    <chr>          <dbl>    <dbl> <dbl>     <dbl>
1 Karitiana.DG USA_Anzick_realigned.SG          Mbuti.DG CHB.DG      0.00379  0.000356 10.6  2.13e- 26
2 Karitiana.DG USA_Anzick_realigned.SG          Mbuti.DG Papuan.DG   0.000825 0.000381  2.17 3.01e-  2
3 Karitiana.DG USA_Anzick_realigned.SG          Mbuti.DG Denisova.DG 0.00981  0.000313 31.3  1.01e-215
4 Karitiana.DG USA_Ancient_Beringian.SG         Mbuti.DG CHB.DG      0.00369  0.000396  9.30 1.35e- 20
5 Karitiana.DG USA_Ancient_Beringian.SG         Mbuti.DG Papuan.DG   0.00112  0.000427  2.63 8.55e-  3
6 Karitiana.DG USA_Ancient_Beringian.SG         Mbuti.DG Denisova.DG 0.00979  0.000341 28.7  1.10e-181
7 Karitiana.DG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG CHB.DG      0.00471  0.000342 13.8  4.74e- 43
8 Karitiana.DG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG Papuan.DG   0.00103  0.000394  2.60 9.23e-  3
9 Karitiana.DG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG Denisova.DG 0.00922  0.000307 30.0  2.58e-198

$rankdrop
# A tibble: 3 × 7
  f4rank   dof    chisq      p dofdiff chisqdiff p_nested
   <int> <int>    <dbl>  <dbl>   <int>     <dbl>    <dbl>
1      2     1    0.463 0.496        3      12.6  0.00553
2      1     4   13.1   0.0109       5    1678.   0      
3      0     9 1691.    0           NA      NA   NA      

```
```
========================================
Resultados qpWave para: Surui.DG 
========================================
$f4_var_rcond
[1] 0.08791833

$f4_var_singular_loadings
NULL

$f4
# A tibble: 9 × 8
  pop1     pop2                             pop3     pop4             est       se     z         p
  <chr>    <chr>                            <chr>    <chr>          <dbl>    <dbl> <dbl>     <dbl>
1 Surui.DG USA_Anzick_realigned.SG          Mbuti.DG CHB.DG      0.00441  0.000356 12.4  2.85e- 35
2 Surui.DG USA_Anzick_realigned.SG          Mbuti.DG Papuan.DG   0.000546 0.000394  1.38 1.66e-  1
3 Surui.DG USA_Anzick_realigned.SG          Mbuti.DG Denisova.DG 0.0105   0.000312 33.7  5.30e-249
4 Surui.DG USA_Ancient_Beringian.SG         Mbuti.DG CHB.DG      0.00431  0.000413 10.4  1.70e- 25
5 Surui.DG USA_Ancient_Beringian.SG         Mbuti.DG Papuan.DG   0.000842 0.000447  1.88 5.96e-  2
6 Surui.DG USA_Ancient_Beringian.SG         Mbuti.DG Denisova.DG 0.0105   0.000339 31.0  1.27e-210
7 Surui.DG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG CHB.DG      0.00533  0.000343 15.5  1.71e- 54
8 Surui.DG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG Papuan.DG   0.000747 0.000385  1.94 5.23e-  2
9 Surui.DG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG Denisova.DG 0.00992  0.000304 32.7  4.35e-234

$rankdrop
# A tibble: 3 × 7
  f4rank   dof    chisq      p dofdiff chisqdiff p_nested
   <int> <int>    <dbl>  <dbl>   <int>     <dbl>    <dbl>
1      2     1    0.447 0.504        3      12.8  0.00505
2      1     4   13.3   0.0101       5    2161.   0      
3      0     9 2174.    0           NA      NA   NA      
```

```
========================================
Resultados qpWave para: PEL.DG 
========================================
$f4_var_rcond
[1] 0.1436458

$f4_var_singular_loadings
NULL

$f4
# A tibble: 9 × 8
  pop1   pop2                             pop3     pop4            est       se     z         p
  <chr>  <chr>                            <chr>    <chr>         <dbl>    <dbl> <dbl>     <dbl>
1 PEL.DG USA_Anzick_realigned.SG          Mbuti.DG CHB.DG      0.00412 0.000314 13.1  3.29e- 39
2 PEL.DG USA_Anzick_realigned.SG          Mbuti.DG Papuan.DG   0.00263 0.000338  7.80 6.30e- 15
3 PEL.DG USA_Anzick_realigned.SG          Mbuti.DG Denisova.DG 0.00871 0.000278 31.4  9.49e-216
4 PEL.DG USA_Ancient_Beringian.SG         Mbuti.DG CHB.DG      0.00402 0.000351 11.4  2.84e- 30
5 PEL.DG USA_Ancient_Beringian.SG         Mbuti.DG Papuan.DG   0.00293 0.000373  7.86 3.93e- 15
6 PEL.DG USA_Ancient_Beringian.SG         Mbuti.DG Denisova.DG 0.00869 0.000297 29.2  5.08e-188
7 PEL.DG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG CHB.DG      0.00504 0.000299 16.8  1.25e- 63
8 PEL.DG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG Papuan.DG   0.00283 0.000333  8.52 1.59e- 17
9 PEL.DG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG Denisova.DG 0.00812 0.000269 30.2  2.72e-200

$rankdrop
# A tibble: 3 × 7
  f4rank   dof    chisq      p dofdiff chisqdiff p_nested
   <int> <int>    <dbl>  <dbl>   <int>     <dbl>    <dbl>
1      2     1    0.548 0.459        3      12.6  0.00571
2      1     4   13.1   0.0108       5    2486.   0      
3      0     9 2499.    0           NA      NA   NA      
```

```
========================================
Resultados qpWave para: MXL.DG 
========================================
$f4_var_rcond
[1] 0.1106277

$f4_var_singular_loadings
NULL

$f4
# A tibble: 9 × 8
  pop1   pop2                             pop3     pop4            est       se     z         p
  <chr>  <chr>                            <chr>    <chr>         <dbl>    <dbl> <dbl>     <dbl>
1 MXL.DG USA_Anzick_realigned.SG          Mbuti.DG CHB.DG      0.00892 0.000332 26.9  1.28e-159
2 MXL.DG USA_Anzick_realigned.SG          Mbuti.DG Papuan.DG   0.00565 0.000362 15.6  5.71e- 55
3 MXL.DG USA_Anzick_realigned.SG          Mbuti.DG Denisova.DG 0.00323 0.000276 11.7  8.11e- 32
4 MXL.DG USA_Ancient_Beringian.SG         Mbuti.DG CHB.DG      0.00883 0.000362 24.4  2.56e-131
5 MXL.DG USA_Ancient_Beringian.SG         Mbuti.DG Papuan.DG   0.00595 0.000379 15.7  1.35e- 55
6 MXL.DG USA_Ancient_Beringian.SG         Mbuti.DG Denisova.DG 0.00322 0.000280 11.5  1.05e- 30
7 MXL.DG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG CHB.DG      0.00985 0.000326 30.2  5.97e-201
8 MXL.DG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG Papuan.DG   0.00585 0.000372 15.7  9.46e- 56
9 MXL.DG USA_Nevada_SpiritCave_11000BP.SG Mbuti.DG Denisova.DG 0.00265 0.000266  9.96 2.29e- 23

$rankdrop
# A tibble: 3 × 7
  f4rank   dof    chisq      p dofdiff chisqdiff p_nested
   <int> <int>    <dbl>  <dbl>   <int>     <dbl>    <dbl>
1      2     1    0.627 0.429        3      10.4   0.0158
2      1     4   11.0   0.0268       5    2175.    0     
3      0     9 2186.    0           NA      NA    NA     
```

Interpretacion:
Al analizar los resultados de qpwave se observa que para las 5 poblaciones target estudiadas (Pima, Karitiana, Surui, PEL, y MXL) se require 3 flujos de ancestría para explicar la relación con las poblaciones de referencia. 


qpadm
**Resultado:**
```
==================================================
 MODELO qpAdm DE 3 VÍAS PARA: Pima.DG 
==================================================

  target  left                             weight    se      z
  <chr>   <chr>                             <dbl> <dbl>  <dbl>
1 Pima.DG USA_Anzick_realigned.SG           -73.3 209.  -0.350
2 Pima.DG USA_Ancient_Beringian.SG           61.3 182.   0.337
3 Pima.DG USA_Nevada_SpiritCave_11000BP.SG   13.0  36.4  0.358


  pat      wt   dof    chisq         p f4rank USA_Anzick_realigned.SG USA_Ancient_Beringian.SG USA_Nevada_SpiritCave_11000BP.SG feasible best  dofdiff chisqdiff p_nested
  <chr> <dbl> <dbl>    <dbl>     <dbl>  <dbl>                   <dbl>                    <dbl>                            <dbl> <lgl>    <lgl>   <dbl>     <dbl>    <dbl>
1 000       0     1    0.490 4.84e-  1      2                   -73.3                     61.3                             13.0 FALSE    NA         NA     NA          NA
2 001       1     2    0.788 6.74e-  1      1                  -151.                     152.                              NA   FALSE    TRUE        0     -7.63        1
3 010       1     2    8.42  1.49e-  2      1                   -27.3                     NA                               28.3 FALSE    TRUE        0     -1.93        1
4 100       1     2   10.3   5.66e-  3      1                    NA                      -35.2                             36.2 FALSE    TRUE       NA     NA          NA
5 011       2     3 1376.    4.02e-298      0                     1                       NA                               NA   TRUE     NA         NA     NA          NA
6 101       2     3 1097.    1.93e-237      0                    NA                        1                               NA   TRUE     NA         NA     NA          NA
7 110       2     3 1341.    1.83e-290      0                    NA                       NA                                1   TRUE     NA         NA     NA          NA
```

```
==================================================
 MODELO qpAdm DE 3 VÍAS PARA: Karitiana.DG 
==================================================

[1] Proporciones de mezcla (weights):
# A tibble: 3 × 5
  target       left                             weight    se      z
  <chr>        <chr>                             <dbl> <dbl>  <dbl>
1 Karitiana.DG USA_Anzick_realigned.SG           -71.5 185.  -0.386
2 Karitiana.DG USA_Ancient_Beringian.SG           59.9 161.   0.372
3 Karitiana.DG USA_Nevada_SpiritCave_11000BP.SG   12.6  33.6  0.374

[2] Ajuste del modelo y subtipos (popdrop):
# A tibble: 7 × 14
  pat      wt   dof    chisq         p f4rank USA_Anzick_realigned.SG USA_Ancient_Beringian.SG USA_Nevada_SpiritCave_11000BP.SG feasible best  dofdiff chisqdiff p_nested
  <chr> <dbl> <dbl>    <dbl>     <dbl>  <dbl>                   <dbl>                    <dbl>                            <dbl> <lgl>    <lgl>   <dbl>     <dbl>    <dbl>
1 000       0     1    0.463 4.96e-  1      2                   -71.5                     59.9                             12.6 FALSE    NA         NA     NA          NA
2 001       1     2    0.763 6.83e-  1      1                  -121.                     122.                              NA   FALSE    TRUE        0     -7.85        1
3 010       1     2    8.61  1.35e-  2      1                   -32.0                     NA                               33.0 FALSE    TRUE        0     -1.99        1
4 100       1     2   10.6   4.98e-  3      1                    NA                      -48.1                             49.1 FALSE    TRUE       NA     NA          NA
5 011       2     3 1473.    4.47e-319      0                     1                       NA                               NA   TRUE     NA         NA     NA          NA
6 101       2     3 1151.    3.32e-249      0                    NA                        1                               NA   TRUE     NA         NA     NA          NA
7 110       2     3 1420.    1.68e-307      0                    NA                       NA                                1   TRUE     NA         NA     NA          NA
```
```
==================================================
 MODELO qpAdm DE 3 VÍAS PARA: Surui.DG 
==================================================

[1] Proporciones de mezcla (weights):
# A tibble: 3 × 5
  target   left                             weight    se      z
  <chr>    <chr>                             <dbl> <dbl>  <dbl>
1 Surui.DG USA_Anzick_realigned.SG           -74.9 184.  -0.407
2 Surui.DG USA_Ancient_Beringian.SG           63.0 160.   0.394
3 Surui.DG USA_Nevada_SpiritCave_11000BP.SG   12.9  34.3  0.377

[2] Ajuste del modelo y subtipos (popdrop):
# A tibble: 7 × 14
  pat      wt   dof    chisq         p f4rank USA_Anzick_realigned.SG USA_Ancient_Beringian.SG USA_Nevada_SpiritCave_11000BP.SG feasible best  dofdiff chisqdiff p_nested
  <chr> <dbl> <dbl>    <dbl>     <dbl>  <dbl>                   <dbl>                    <dbl>                            <dbl> <lgl>    <lgl>   <dbl>     <dbl>    <dbl>
1 000       0     1    0.447 5.04e-  1      2                   -74.9                     63.0                             12.9 FALSE    NA         NA     NA          NA
2 001       1     2    0.745 6.89e-  1      1                  -113.                     114.                              NA   FALSE    TRUE        0     -8.01        1
3 010       1     2    8.76  1.26e-  2      1                   -38.8                     NA                               39.8 FALSE    TRUE        0     -1.94        1
4 100       1     2   10.7   4.75e-  3      1                    NA                      -68.6                             69.6 FALSE    TRUE       NA     NA          NA
5 011       2     3 1809.    0              0                     1                       NA                               NA   TRUE     NA         NA     NA          NA
6 101       2     3 1351.    1.17e-292      0                    NA                        1                               NA   TRUE     NA         NA     NA          NA
7 110       2     3 1709.    0              0                    NA                       NA                                1   TRUE     NA         NA     NA          NA
```

```
==================================================
 MODELO qpAdm DE 3 VÍAS PARA: PEL.DG 
==================================================

[1] Proporciones de mezcla (weights):
# A tibble: 3 × 5
  target left                             weight    se      z
  <chr>  <chr>                             <dbl> <dbl>  <dbl>
1 PEL.DG USA_Anzick_realigned.SG           -89.8 352.  -0.255
2 PEL.DG USA_Ancient_Beringian.SG           76.1 307.   0.248
3 PEL.DG USA_Nevada_SpiritCave_11000BP.SG   14.7  55.1  0.267

[2] Ajuste del modelo y subtipos (popdrop):
# A tibble: 7 × 14
  pat      wt   dof    chisq         p f4rank USA_Anzick_realigned.SG USA_Ancient_Beringian.SG USA_Nevada_SpiritCave_11000BP.SG feasible best  dofdiff chisqdiff p_nested
  <chr> <dbl> <dbl>    <dbl>     <dbl>  <dbl>                   <dbl>                    <dbl>                            <dbl> <lgl>    <lgl>   <dbl>     <dbl>    <dbl>
1 000       0     1    0.548 4.59e-  1      2                   -89.8                     76.1                             14.7 FALSE    NA         NA     NA          NA
2 001       1     2    0.806 6.68e-  1      1                  -213.                     214.                              NA   FALSE    TRUE        0     -7.94        1
3 010       1     2    8.74  1.26e-  2      1                   -30.8                     NA                               31.8 FALSE    TRUE        0     -1.71        1
4 100       1     2   10.5   5.36e-  3      1                    NA                      -37.4                             38.4 FALSE    TRUE       NA     NA          NA
5 011       2     3 1236.    1.17e-267      0                     1                       NA                               NA   TRUE     NA         NA     NA          NA
6 101       2     3 1032.    1.89e-223      0                    NA                        1                               NA   TRUE     NA         NA     NA          NA
7 110       2     3 1294.    3.67e-280      0                    NA                       NA                                1   TRUE     NA         NA     NA          NA
```

```
==================================================
 MODELO qpAdm DE 3 VÍAS PARA: MXL.DG 
==================================================

[1] Proporciones de mezcla (weights):
# A tibble: 3 × 5
  target left                             weight      se        z
  <chr>  <chr>                             <dbl>   <dbl>    <dbl>
1 MXL.DG USA_Anzick_realigned.SG           504.  201528.  0.00250
2 MXL.DG USA_Ancient_Beringian.SG         -433.  174096. -0.00249
3 MXL.DG USA_Nevada_SpiritCave_11000BP.SG  -70.5  27444. -0.00257

[2] Ajuste del modelo y subtipos (popdrop):
# A tibble: 7 × 14
  pat      wt   dof    chisq         p f4rank USA_Anzick_realigned.SG USA_Ancient_Beringian.SG USA_Nevada_SpiritCave_11000BP.SG feasible best  dofdiff chisqdiff p_nested
  <chr> <dbl> <dbl>    <dbl>     <dbl>  <dbl>                   <dbl>                    <dbl>                            <dbl> <lgl>    <lgl>   <dbl>     <dbl>    <dbl>
1 000       0     1    0.627 4.29e-  1      2                   504.                    -433.                             -70.5 FALSE    NA         NA     NA          NA
2 001       1     2    0.802 6.70e-  1      1                  -132.                     133.                              NA   FALSE    TRUE        0     -6.29        1
3 010       1     2    7.09  2.88e-  2      1                    15.7                     NA                              -14.7 FALSE    TRUE        0     -1.50        1
4 100       1     2    8.59  1.36e-  2      1                    NA                       14.4                            -13.4 FALSE    TRUE       NA     NA          NA
5 011       2     3  882.    7.59e-191      0                     1                       NA                               NA   TRUE     NA         NA     NA          NA
6 101       2     3  725.    9.80e-157      0                    NA                        1                               NA   TRUE     NA         NA     NA          NA
7 110       2     3 1078.    1.90e-233      0                    NA                       NA                                1   TRUE     NA         NA     NA          NA
```


Finalmente al analizar el resultado de qpAdmixture se evidencia lo siguiente. Destacan en todas las poblaciones target estudiadas errores estandar muy elevados, lo que tmabien en ciertos casos se correlacona con adm_weignts sumammente altos. SIn duda este resultados llama la atención y orienta a una posibilidad de cambiar las poblaciones de referencia o outgroups para justamente lograr determinar las fuentes de mestizaje de las poblaciones latinoamericanas modernas estudiadas.


## Proyecto 4: Medieval / Iberian admixture

No se pudo trabajar con la población Onge, y para escoger población del norte de áfrica se designo la población Marrueca y Tunes_africa. 

Pregunta de investigación: ¿Tiene la poblacion mediterrania-ibérica mestizaje con población del norte de áfrica?.

Objetivo de estudio: Evaluar si la población IBS posee mestizaje con poblaciones Iberomaurusian e YRI, las cuales corresponden a poblaciones ancestrales del norte de áfrica. Como población outgroup se utilizó la población Mbuti. 

Para el calculo de f2, principalmetne por limitaciones de RAM del computador, solo se calculo para un numero limitado de poblacioens (IBS, Ibermaurusian, YRI y Mbuti).

Calculo F3:
**Resultado:**
```
  pop1     pop2   pop3                         est       se     z     p
  <chr>    <chr>  <chr>                      <dbl>    <dbl> <dbl> <dbl>
1 Mbuti.DG IBS.DG Morocco_Iberomaurusian.AG 0.0547 0.000666  82.2     0
2 Mbuti.DG IBS.DG YRI.DG                    0.0253 0.000392  64.5     0
```

De estos resultados se desprende una mayor proporción de deriva génica entre la población IBS y Morocoo_Iberomaurusian. De igual forma un porcentaje de la estructura génomica de IBS también contiene alelos comaportidos con la población YRI.DG, provieniente esta última de Ibadán, Nigeria. 

Calculo F4:
**Resultado:**
```
  pop1   pop2                      pop3                      pop4          est       se       z        p
  <chr>  <chr>                     <chr>                     <chr>       <dbl>    <dbl>   <dbl>    <dbl>
1 IBS.DG Morocco_Iberomaurusian.AG Morocco_Iberomaurusian.AG Mbuti.DG -0.0518  0.000442 -117.   0       
2 IBS.DG Morocco_Iberomaurusian.AG YRI.DG                    Mbuti.DG -0.00126 0.000140   -8.98 2.71e-19
3 IBS.DG YRI.DG                    Morocco_Iberomaurusian.AG Mbuti.DG  0.0282  0.000397   71.0  0       
4 IBS.DG YRI.DG                    YRI.DG                    Mbuti.DG -0.00193 0.000128  -15.1  9.94e-52
```

Este resultado indica que la población ibérica IBS.DG comaprte una mayor cantidad de alelos con la población Morocco_Iberomaurisian en comapración a la poblción Yri. 


Calculo qpwave:

**Resultado:**
```
========================================
Resultados qpWave para: IBS.DG 
========================================
$f4_var_rcond
[1] 0.3989661

$f4_var_singular_loadings
NULL

$f4
# A tibble: 2 × 8
  pop1   pop2                      pop3     pop4          est       se     z         p
  <chr>  <chr>                     <chr>    <chr>       <dbl>    <dbl> <dbl>     <dbl>
1 IBS.DG Morocco_Iberomaurusian.AG Mbuti.DG Papuan.DG -0.0110 0.000380 -28.8 1.22e-182
2 IBS.DG YRI.DG                    Mbuti.DG Papuan.DG -0.0347 0.000420 -82.5 0        

$rankdrop
# A tibble: 2 × 7
  f4rank   dof       chisq     p dofdiff chisqdiff p_nested
   <int> <int>       <dbl> <dbl>   <int>     <dbl>    <dbl>
1      1     0    0.000105    NA       2     6808.        0
2      0     2 6808.           0      NA       NA        NA
```

De estos resultados se puede concluir que se dependen de al menos de >2 flujos de ancestría para explciar la relación entre la población target y las poblacoensde referencia. Vale destacar que el numero de poblaciones que se pudieron incluir fueron limitadas por la capacidad del computador de procesar f2, por lo que los resultados presentan dicho sesgo. 

Calculo qpadmin

**Resultado:**
```
==================================================
 MODELO qpAdm PARA TARGET: IBS.DG 
==================================================

[1] Proporciones de mezcla (weights):
# A tibble: 2 × 5
  target left                      weight     se     z
  <chr>  <chr>                      <dbl>  <dbl> <dbl>
1 IBS.DG Morocco_Iberomaurusian.AG  1.46  0.0240  60.9
2 IBS.DG YRI.DG                    -0.462 0.0240 -19.2

[2] Evaluación de ajuste global y submodelos (popdrop):
# A tibble: 3 × 13
  pat      wt   dof       chisq          p f4rank Morocco_Iberomaurusian.AG YRI.DG feasible best  dofdiff chisqdiff p_nested
  <chr> <dbl> <dbl>       <dbl>      <dbl>  <dbl>                     <dbl>  <dbl> <lgl>    <lgl>   <dbl>     <dbl>    <dbl>
1 00        0     0    0.000105 NA              1                      1.46 -0.462 FALSE    NA         NA       NA        NA
2 01        1     1  960.        9.63e-211      0                      1    NA     TRUE     TRUE        0    -6907.        1
3 10        1     1 7866.        0              0                     NA     1     TRUE     TRUE       NA       NA        NA
```

AL evidenciar valores de weight_adm se observan valores poco confiables. Esto tambien se puede deber al kimitado grupo de poblaciones de referecnia que se incluyeron en el estudio. No es posible desprender una conclusión acerca del linaje de la pobalción IBS. ES posible que su genoma sea comuesto también por un linaje de otras poblacioens no incluidas en este estudio. 



---

**Conclusiones proyecto:**
 
En el presente trabajo se praćtico el uso de estadisticos f, calculo de qpwave y qpadmin, con el fin final de evaluar preguntas científicas acorde a la genetica de poblaciones, y particularmente, rsponder a como se compone el genoma de las distitnas poblaciones modernas y ancestrales, como ocurren los procesos de ancestria y mestizaje, y finalemnte como se van formando las poblaciones a lo largo del tiempo. A modo personal me parecío un tema sumamente interesante, permite responder con evidencia robusta interrogantes naturales que ocurren en poblaciones. La idea de saber como fue formado y de dondo proviene la ancestria o nos permite reconstruir la historia del poblamiento humano y sus complejos eventos de mestizaje, orientando a su vez la comprensión de la diversidad genética actual.
Como principales desafios de esta tarea principalmente fue entender la interpretación de los diversos test. Destaco comandos de admintools faciles de seguir por lo que sta vez el trabajo bioinformático no fue la gran limitante, sino entender justamente que significaba, por ejemplo, valores de f3 o f4 (+) o (-). Otra gran limitante ya se relacionó con la memoria RAM del comutador. Lamentablemente no fue posible realizar en su totalidad los calculos de f2 para la totalidad de poblaciones que me hubiera gustado incluir para responder las preguntas, lo que se instala como una de las pirncipales limitantes que presentó el desarrollo de los distintos proyectos. 


Referencias:
1. Peter B. M. (2016). Admixture, Population Structure, and F-Statistics. Genetics, 202(4), 1485–1501. https://doi.org/10.1534/genetics.115.183913
2. Mallick, S., Micco, A., Mah, M., Ringbauer, H., Lazaridis, I., Olalde, I., Patterson, N., & Reich, D. (2024). The Allen Ancient DNA Resource (AADR) a curated compendium of ancient human genomes. Scientific data, 11(1), 182. https://doi.org/10.1038/s41597-024-03031-7
