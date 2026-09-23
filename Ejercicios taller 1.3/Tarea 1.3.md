# Tarea sesión 3
# **_Bioinformatica y analisis genómicos reproducibles_**
Tarea Sesión 1.3
- Autor: Joaquín Ignacio Torres Núñez
- Curso: Bioinformatics and Reproducible Genomic Analysis
- Profesor: Ricardo Verdugo
- Unidad 1 – Sesion 3

## Ejercicio 1 crea una variable con el logaritmo base 10 de 50 y súmalo a otra variable cuyo valor sea igual a 5.
````
log
x=log10(50)
y=5
z=x + y
print(z)
````
## Ejercicio 2 suma el número 2 a todos los números entre 1 y 150.
```
x=(1:150)
y=x+2
print(y)
```

## Ejercicio 3 ¿cuántos números son mayores a 20 en el vector -13432:234?
```
?sum
vector=c(-13432:234)
sum(vector>20)
```

## Ejercicio 4 Carga en R el archivo PracUni1Ses3/maices/meta/maizteocintle_SNP50k_meta_extended.txt y ponlo en un objeto de R llamado meta_maiz.
```
getwd()
?read.delim
# Cargar el archivo de metadatos
meta_maiz <- read.delim("./BioinfinvRepro-master/BioinfinvRepro-master/Unidad1/Sesion3/PracUni1Ses3/maices/meta/maizteocintle_SNP50k_meta_extended.txt")
````

## Ejercico 5 
### ejemplo
```
for (i in 2:10){
  print(paste(i, "elefantes se columpiaban sobre la tela de una araña"))
}
```
### Ejercicio 5.1 Escribe un for loop para que divida 35 entre 1:10 e imprima el resultado en la consola.
```
for (i in 1:10) {
  x = 35/i
  print(paste("El resultado de la división 35 /", i, "es igual a", x))
}
```

### Ejercicio 5.2 Modifica el loop anterior para que haga las divisiones solo para los números nones (con un comando, NO con c(1,3,...)). Pista: next.
```
for (i in 1:10) {
  if (i %% 2 == 0) {
    next
  }
  x <- 35/i
  print(paste("El resultado de la división 35 /", i, "es igual a", x))
}
```

### Ejercicio 5.3 Modifica el loop anterior para que los resultados de correr todo el loop se guarden en una df de dos columnas, la primera debe tener el texto "resultado para x" (donde x es cada uno de los elementos del loop) y la segunda el resultado correspondiente a cada elemento del loop. Pista: el primer paso es crear un vector fuera del loop. Ejemplo:
```
# 1. Crear un data frame vacío fuera del loop
df_resultados <- data.frame()

# 2. Ejercicio del loop (usando 'next' para filtrar solo números nones)
for (i in 1:10) {
  x <- 35 / i
  # 3. Ir agregando cada nueva fila con rbind()
  df_resultados <- rbind(df_resultados, data.frame(
    columna1 = paste("resultado para", i),
    columna2 = x
  ))
}
# 4. Mostrar la matriz resultante
df_resultados
```
### Ejercicio 6 Abre en RStudio el script PracUni1Ses3/mantel/bin/1.IBR_testing.r. Este script realiza un análisis de aislamiento por resistencia con Fst calculadas con ddRAD en Berberis alpina. Lee el código del script y determina:

#### ¿qué hacen los dos for loops del script?: 
- El primer loop realiza un ciclo for a traves de los 13 vectores calculando para cada vector la matriz de distancia efectiva.
````
 for(i in c("present", "ccsm", "miroc", "flat", "1800", "2000", "2300", "2500", "2700", "3000", "3300", "3500", "4000")) {
    
    ## define resistances.out files
    resfile <- paste0(circfolder, "/Balpina_", i, "_resistances.out")
      
    ### Get effective distances
    
    eff.dist<-read.effdist(file=resfile, popNames=popNamesFP, des.order=popNames)

````
- El segundo loop itera por cada vector calculando el test de mantel, donde compara la matriz de distancia ambiental (present, csm, miroc, etc) con la matriz genetica (B.FstLin)

````
 for(i in c("present", "ccsm", "miroc", "flat", "1800", "2000", "2300", "2500", "2700", "3000", "3300", "3500", "4000")) {
  
    print(paste("Results for", i))
    
    # Mantel test 
    print("Mantel test")
    x<-mantel.rtest(as.dist(get(paste0("B.",i))), as.dist(B.FstLin), nrepet=10000)
    print(x)

````


#### ¿qué paquetes necesitas para correr el script?: 
- (ade4), (ggplop), (sp)
#### ¿qué archivos necesitas para correr el script?: 
- (read.fst_summary_fix.R), (read.effdist.R), (surveyed_mountains.tsv), (BerSS.sumstats.tsv), (Balpina_focalpoints.txt), 

### Ejercicio 7 Escribe una función llamada calc.tetha que te permita calcular tetha dados Ne y u como argumentos. Recuerda que tetha =4Neu.
```
#Red Bananas
calc.tetha=function(Ne, u){
  tetha=4*Ne*u
  return(tetha)
}
```
### Ejercicio 8 Al script del ejercicio de las pruebas de Mantel, agrega el código necesario para realizar un Partial Mantel test entre la matriz Fst, y las matrices del presente y el LGM, parcializando la matriz flat. Necesitarás el paquete vegan. Include the comment: Elefante blanco

```
## Alicia Mastretta Yanes
## Perform Mantel tests with different resistance surfaces to test for IBR in B. alpina
#Partial mantel test permite el analisi de 3 matrices, a diferencia de mantel test. estima la correlacion de 2 matrices (A y B), mientras se controla el efecto de la matriz C. Se suele utilizar como matriz C una atriz de distancia derivado de un parametro ambiental. 

# start with a fresh brain
rm(list = ls())

# Load libraries, carga las librerias necesarias, se incluye libreria vegan, requerido para Partial mantel
library(ade4)
library(ggplot2)
library(sp)
library(vegan)

########### Get data  ###############
## Source home made funcs to load data, carga funciones necesarias
source("read.fst_summary_fix.R")
source("read.effdist.R")

### Geographic ###

# For reference, get Population ID codes as used for Circuitscape and Map plotting
# this is not the same order than PopuplationMaps used for Stacks. Careful. 
points.info<-read.delim("../spatial/surveyed_mountains.tsv")
points.info
points.infoxy<-as.matrix(points.info[,c(5,6)], longlat=TRUE)

# Get geographic distances. Carga coordenadas geograficas
GeoDist<-spDists(points.infoxy, ) 
colnames(GeoDist)<-points.info$Key
rownames(GeoDist)<-points.info$Key

### define paths for loading data
genfolder<-"../genetic"
circfolder<-"../spatial/resdist"

### Genetic ###

## define pop names as in Stack population maps 
# check pop map order
readLines(paste0(genfolder, "/BerSS.sumstats.tsv"), n=6)
# define popNames
popNames=c("Aj","Iz","Ma","Pe","Tl","To") 

## Get Fst pairwise matrix. Se guarda variable B.Fst que contiene una matriz de distancia genetica entre los popNames
  B.Fst <-read.fst_summary_fix(file=paste0(genfolder, "/BerSS.fst_summary.tsv"),
                         popNames=popNames)   


### Effective distances ###

### Get general info and paths
  ## define pop names as in Circuitscape focal points order 
  # get focal points. Verificva si esta de forma correcta el orden dre los puntos focales de las muestras
  focpoints<-read.delim(paste0(circfolder, "/Balpina_focalpoints.txt"), header=FALSE)
  # get info of focal points
  x<-points.info$ID %in% focpoints[,1] 
  focpoints<- points.info[x,] 
  # get PopNames in order of focal points
  popNamesFP<-as.vector(focpoints$Key)
  popNamesFP
  # get PopNames in order of Stacks output PopMap
  popNames=popNames
  popNames
 

### Get effective distance matrix and mean of it for each raster 

  for(i in c("present", "ccsm", "miroc", "flat", "1800", "2000", "2300", "2500", "2700", "3000", "3300", "3500", "4000")) {
    
    ## define resistances.out files
    resfile <- paste0(circfolder, "/Balpina_", i, "_resistances.out")
      
    ### Get effective distances
    
    eff.dist<-read.effdist(file=resfile, popNames=popNamesFP, des.order=popNames)
    
    ### Estimate mean effective distance by population
    mean.effD <- apply(eff.dist, 2, mean)  
      
    ### Name output data
    assign(paste0("B.", i), eff.dist)  # effective distance mat
    assign(paste0("B.mean.", i), mean.effD) # mean effective distances
  
  }    
    
### Get Geographic distances for this spp localities
B.GeoDist<-GeoDist[match(popNames,rownames(GeoDist)), match(popNames,colnames(GeoDist))] #get right order
B.GeoDist


###########  Isolation by Resistance ###############

##### Perform Mantel test between the Fst matrix and the present and LGM effective distances
# function for lm and plotting
source("DistPlot.R")


### Berberis
  # Linearize as suggested by Rousset (1997) for IBD using FST/(1 - FST)
  B.FstLin<- B.Fst/(1-B.Fst)
  
  # run mantel test for each condition
IBRresults<-c("rster", "MTpvalue", "MTr")
  for(i in c("present", "ccsm", "miroc", "flat", "1800", "2000", "2300", "2500", "2700", "3000", "3300", "3500", "4000")) {
  
    print(paste("Results for", i))
    
    # Mantel test 
    print("Mantel test")
    x<-mantel.rtest(as.dist(get(paste0("B.",i))), as.dist(B.FstLin), nrepet=10000)
    print(x)
    
    # Plot
    DistPlot(get(paste0("B.",i)), B.FstLin, plotnames=FALSE,
            ylabel=expression("F"[ST]*"/(1 - "[FST]*")"), xlabel=paste("Effective distance", i))
               
    # get info for df  
    MTpvalue<-round(x$pvalue, 6)
    MTr<-round(x$obs, 4)
    
    # put results in dataframe
    rster<-paste(i)
    IBRresults<-rbind(IBRresults, c(rster, MTpvalue, MTr))
  } 

IBRresults

#### PARTIAL MANTEL TEST #######
# Elefante blanco

fst_d  <- as.dist(B.FstLin)  # Matriz genética
flat_d <- as.dist(B.flat)    # Matriz control/distancia geografica

# 1. Modelo Presente
res_present <- mantel.partial(fst_d, as.dist(B.present), flat_d, method = "pearson", permutations = 10000)
print(res_present)

# 2. Modelo LGM CCSM
res_ccsm <- mantel.partial(fst_d, as.dist(B.ccsm), flat_d, method = "pearson", permutations = 10000)
print(res_ccsm)

# 3. Modelo LGM MIROC
res_miroc <- mantel.partial(fst_d, as.dist(B.miroc), flat_d, method = "pearson", permutations = 10000)
print(res_miroc)


## session info
sessionInfo()

```

### Ejercicio 9 Escribe un script que debe estar guardado en PracUni1Ses3/maices/bin y llamarse ExplorandoMaiz.R, que 1) cargue en R el archivo PPracUni1Ses3maices/meta/maizteocintle_SNP50k_meta_extended.txt y 2) responda lo siguiente.

## ==============================================================================
## Script: ExplorandoMaiz.R
## Ubicación: PracUni1Ses3/maices/bin/ExplorandoMaiz.R
## Autor: Joaquín Ignacio Torres Núñez
## Descripción: Carga y exploración de metadatos de maíz y teocintle
## ==============================================================================
```
library(dplyr)
maiz.data=read.delim("./BioinfinvRepro-master/BioinfinvRepro-master/Unidad1/Sesion3/PracUni1Ses3/maices/meta/maizteocintle_SNP50k_meta_extended.txt")
```
- ¿Qué tipo de objeto creamos al cargar la base?
#funcion class entrega el tipo de archivo creado
```
class(maiz.data)
```
- ¿Cómo se ven las primeras 6 líneas del archivo?
#función head entrega por defecto las primeras 6 filas, si se quisiera seleccionar mas filas se agrega head(...,numero de filas)
```
head(maiz.data,6)
```
- ¿Cuántas muestras hay?
#paquete de dplyr, entrega en numero total de filas
```
count(maiz.data)
```
- ¿De cuántos estados se tienen muestras?
#data: La tabla de datos (data frame) con la que vas a trabajar.sintaxis es (.data,), se agrega la columna o lista de columnas por las cuales quieres agrupar y contar.
```
count(maiz.data, Estado)
nrow(count(maiz.data, Estado))
```
- ¿Cuántas muestras fueron colectadas antes de 1980?
#muestra el numero de filas >1980 en su año de colecta, crea un data frame, TRUE es el valor que se busca)
```
count(maiz.data,A.o._de_colecta<1980)
```
- ¿Cuántas muestras hay de cada raza?
#Se indica cuentas muestras hay de cada raza
```
count(maiz.data,Raza)
```
- En promedio ¿a qué altitud fueron colectadas las muestras?
#Promedio de alturas, summarise() permite realizar calculos arimeticos a columnas definidas
```
count(maiz.data,Altitud)
?summarise
summarise(maiz.data,mean(Altitud, na.rm = TRUE))
```
- ¿Y a qué altitud máxima y mínima fueron colectadas?
```
#maximos y minimos
summarise(maiz.data,paste('La altitud máxima es:',max(Altitud, na.rm = TRUE)))  
summarise(maiz.data,paste('La altitud mínima es:',min(Altitud, na.rm = TRUE)))
```
- Crea una nueva df de datos sólo con las muestras de la raza Olotillo
#nuevo data frame
```
df_olotillo = filter(maiz.data, Raza == "Olotillo")
```

- Crea una nueva df de datos sólo con las muestras de la raza Reventador, Jala y Ancho
#nuevo data frame de mas de 1 raza
```
df_mix = filter(maiz.data, Raza == "Reventador"| Raza == "Jala"|Raza == "Ancho")
dim(df_mix)
```
- Escribe la matriz anterior a un archivo llamado "submat.cvs" en /meta

#Necesario importar libreria readr
```
library(readr)

write_csv(df_mix, "./BioinfinvRepro-master/BioinfinvRepro-master/Unidad1/Sesion3/PracUni1Ses3/maices/meta/submatJT.csv")
```
### Ejercicio 10 
El resultado del Ejercicio 10 se encuentra en [Wiki sesion 1.3] (https://github.com/joaquintorresnunez/Tareas_BioinfRepro2026_JITN/wiki/Notas-sesion-1.3)


## Conclusion




