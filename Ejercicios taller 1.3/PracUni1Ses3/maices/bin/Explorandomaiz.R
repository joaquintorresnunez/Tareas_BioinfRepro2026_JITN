### Ejercicio 9 Escribe un script que debe estar guardado en PracUni1Ses3/maices/bin y llamarse ExplorandoMaiz.R, que 1) cargue en R el archivo PPracUni1Ses3maices/meta/maizteocintle_SNP50k_meta_extended.txt y 2) responda lo siguiente.

## ==============================================================================
## Script: ExplorandoMaiz.R
## Ubicación: PracUni1Ses3/maices/bin/ExplorandoMaiz.R
## Autor: Joaquín Ignacio Torres Núñez
## Descripción: Carga y exploración de metadatos de maíz y teocintle
## ==============================================================================

##### 1. Carga de librerias#####

library(dplyr)
library(readr)

##### 2. Carga de datos#####

maiz.data=read.delim("./BioinfinvRepro-master/Unidad1/Sesion3/PracUni1Ses3/maices/meta/maizteocintle_SNP50k_meta_extended.txt")

##### 3. Evaluación formato y contenido del archivo##### 

#funcion class entrega el tipo de archivo creado
class(maiz.data)

#función head entrega por defecto las primeras 6 filas, si se quisiera seleccionar mas filas se agrega head(...,numero de filas)
head(maiz.data,6)

#paquete de dplyr, entrega en numero total de filas, lo que se traduce en el numero de muestras

count(maiz.data)

##### 4. Conteo y agrupaciones por variables#####

#Conteo de muestras por Estado
count(maiz.data, Estado)
nrow(count(maiz.data, Estado))

#Conteo de muestras colectadas antes del año 1980 (retorna TRUE/FALSE)

count(maiz.data,A.o._de_colecta<1980)

#Conteo de cuántas muestras hay por cada Raza
count(maiz.data,Raza)

####Análisis Estadístico de Altitud####

#Promedio de alturas, summarise() permite realizar calculos arimeticos a columnas definidas
count(maiz.data,Altitud)
?summarise
summarise(maiz.data,mean(Altitud, na.rm = TRUE))

#maximos y minimos
summarise(maiz.data,paste('La altitud máxima es:',max(Altitud, na.rm = TRUE)))  
summarise(maiz.data,paste('La altitud mínima es:',min(Altitud, na.rm = TRUE)))

####Filtrado y Creación de Subconjuntos (Data Frames)####

#nuevo data frame con las muestras de la raza Olotillo

df_olotillo = filter(maiz.data, Raza == "Olotillo")

#nuevo data frame de mas de 1 raza Reventador, Jala y Ancho

df_mix = filter(maiz.data, Raza == "Reventador"| Raza == "Jala"|Raza == "Ancho")
dim(df_mix)

####Exportación de Resultados####

#Escribir la matriz anterior a un archivo llamado "submat.cvs" en /meta, Necesario importar libreria readr

write_csv(df_mix, "./BioinfinvRepro-master/BioinfinvRepro-master/Unidad1/Sesion3/PracUni1Ses3/maices/meta/submatJT.csv")