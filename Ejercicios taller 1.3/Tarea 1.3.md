# Tarea sesión 3
# **_Bioinformatica y analisis genómicos reproducibles_**
Tarea Sesión 1.3
- Autor: Joaquín Ignacio Torres Núñez
- Curso: Bioinformatics and Reproducible Genomic Analysis
- Profesor: Ricardo Verdugo
- Unidad 1 – Sesion 3

## Objetivo general:
- Desarrollar habilidades fundamentales de programación en R y aplicar metodologías bioinformáticas para la manipulación y análisis estadístico de datos.

## Objetivos específicos:
- Implementar ciclos iterativos y funciones personalizadas en R.
- Manipular metadata con paquetes en R
- Aplicar anotaciones de códigos y manejo de sripts en R.

## Ejercicio 1: Crea una variable con el logaritmo base 10 de 50 y súmalo a otra variable cuyo valor sea igual a 5.
````
log
x=log10(50)
y=5
z=x + y
print(z)
````
## Ejercicio 2: Suma el número 2 a todos los números entre 1 y 150.
```
x=(1:150)
y=x+2
print(y)
```

## Ejercicio 3: ¿Cuántos números son mayores a 20 en el vector -13432:234?
```
?sum
vector=c(-13432:234)
sum(vector>20)
```

## Ejercicio 4: Carga en R el archivo PracUni1Ses3/maices/meta/maizteocintle_SNP50k_meta_extended.txt y ponlo en un objeto de R llamado meta_maiz.
```
getwd()
?read.delim
# Cargar el archivo de metadatos
meta_maiz <- read.delim("./BioinfinvRepro-master/BioinfinvRepro-master/Unidad1/Sesion3/PracUni1Ses3/maices/meta/maizteocintle_SNP50k_meta_extended.txt")
````

## Ejercicio 5: 
### Ejercicio 5.1: Escribe un for loop para que divida 35 entre 1:10 e imprima el resultado en la consola.
```
for (i in 1:10) {
  x = 35/i
  print(paste("El resultado de la división 35 /", i, "es igual a", x))
}
```

### Ejercicio 5.2: Modifica el loop anterior para que haga las divisiones solo para los números nones (con un comando, NO con c(1,3,...)). Pista: next.
```
for (i in 1:10) {
  if (i %% 2 == 0) {
    next
  }
  x <- 35/i
  print(paste("El resultado de la división 35 /", i, "es igual a", x))
}
```

### Ejercicio 5.3: Modifica el loop anterior para que los resultados de correr todo el loop se guarden en una df de dos columnas, la primera debe tener el texto "resultado para x" (donde x es cada uno de los elementos del loop) y la segunda el resultado correspondiente a cada elemento del loop. Pista: el primer paso es crear un vector fuera del loop. Ejemplo:
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
### Ejercicio 6: Abre en RStudio el script PracUni1Ses3/mantel/bin/1.IBR_testing.r. Este script realiza un análisis de aislamiento por resistencia con Fst calculadas con ddRAD en Berberis alpina. Lee el código del script y determina:

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

### Ejercicio 7: Escribe una función llamada calc.tetha que te permita calcular tetha dados Ne y u como argumentos. Recuerda que tetha =4Neu.
```
#Red Bananas
calc.tetha=function(Ne, u){
  tetha=4*Ne*u
  return(tetha)
}
```
### Ejercicio 8: Al script del ejercicio de las pruebas de Mantel, agrega el código necesario para realizar un Partial Mantel test entre la matriz Fst, y las matrices del presente y el LGM, parcializando la matriz flat. Necesitarás el paquete vegan. Include the comment: Elefante blanco

```
## Alicia Mastretta Yanes
## Perform Mantel tests with different resistance surfaces to test for IBR in B. alpina
#Partial mantel test permite el analisi de 3 matrices, a diferencia de mantel test. estima la correlacion de 2 matrices (A y B), mientras se controla el efecto de la matriz C. Se suele utilizar como matriz C una atriz de distancia derivado de un parametro ambiental. 

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

### Ejercicio 9: Escribe un script que debe estar guardado en PracUni1Ses3/maices/bin y llamarse ExplorandoMaiz.R, que 1) cargue en R el archivo PPracUni1Ses3maices/meta/maizteocintle_SNP50k_meta_extended.txt y 2) responda lo siguiente.

1. Script ExplorandoMaiz.R disponible en [Ejercicios taller 1.3/PracUni1Ses3/maices/bin](https://github.com/joaquintorresnunez/Tareas_BioinfRepro2026_JITN/blob/0d0fe3a3bc3848bd2bcf524302b7f75dfa6ad2df/Ejercicios%20taller%201.3/PracUni1Ses3/maices/bin/Explorandomaiz.R)
```
library(dplyr)
maiz.data=read.delim("./BioinfinvRepro-master/BioinfinvRepro-master/Unidad1/Sesion3/PracUni1Ses3/maices/meta/maizteocintle_SNP50k_meta_extended.txt")
```
2. Responda lo siguiente:  
- ¿Qué tipo de objeto creamos al cargar la base?: Entrega un data.frame
- ¿Cómo se ven las primeras 6 líneas del archivo?: función head entrega por defecto las primeras 6 filas, si se quisiera seleccionar mas filas se agrega head(...,numero de filas)
```
head(maiz.data,6)
```
- ¿Cuántas muestras hay?: paquete de dplyr, entrega en numero total de filas, son en total 165 muestras
- ¿De cuántos estados se tienen muestras?: 19 estados.
- ¿Cuántas muestras fueron colectadas antes de 1980?: 8 muestars recolectadas antes de 1980.
- ¿Cuántas muestras hay de cada raza?:
  | # | Raza | n |
| :---: | :--- | ---: |
| 1 | Ancho | 3 |
| 2 | Apachito | 2 |
| 3 | Arrocillo | 4 |
| 4 | Azul | 2 |
| 5 | Blando de Sonora | 1 |
| 6 | Bofo | 1 |
| 7 | Cacahuacintle | 5 |
| 8 | Celaya | 3 |
| 9 | Chalqueño | 7 |
| 10 | Chapalote | 2 |
| 11 | Comiteco | 5 |
| 12 | Complejo Serrano de Jalisco | 2 |
| 13 | Conejo | 4 |
| 14 | Coscomatepec | 3 |
| 15 | Cristalino de Chihuahua | 2 |
| 16 | Cónico | 16 |
| 17 | Cónico Norteño | 3 |
| 18 | Dulce | 1 |
| 19 | Dulcillo del Noroeste | 2 |
| 20 | Dzit-Bacal | 3 |
| 21 | Elotero de Sinaloa | 5 |
| 22 | Elotes Cónicos | 14 |
| 23 | Elotes Occidentales | 4 |
| 24 | Gordo | 2 |
| 25 | Jala | 4 |
| 26 | Mushito | 3 |
| 27 | Nal-tel de Altura | 5 |
| 28 | Olotillo | 6 |
| 29 | Olotón | 4 |
| 30 | Onaveño | 2 |
| 31 | Palomero Toluqueño | 1 |
| 32 | Palomero de Chihuahua | 1 |
| 33 | Pepitilla | 4 |
| 34 | Ratón | 3 |
| 35 | Reventador | 2 |
| 36 | Tablilla de Ocho | 2 |
| 37 | Tabloncillo | 4 |
| 38 | Tabloncillo Perla | 3 |
| 39 | Tehua | 2 |
| 40 | Tepecintle | 4 |
| 41 | Tuxpeño | 4 |
| 42 | Tuxpeño Norteño | 2 |
| 43 | Vandeño | 4 |
| 44 | Zamorano Amarillo | 3 |
| 45 | Zapalote Chico | 1 |
| 46 | Zapalote Grande | 1 |
| 47 | Zea m. mexicana | 2 |
| 48 | Zea m. parviglumis | 2 |

- En promedio ¿a qué altitud fueron colectadas las muestras?: altitud media  1519.242
- ¿Y a qué altitud máxima y mínima fueron colectadas?: La altitud máxima es: 2769,  La altitud mínima es: 5
z.data,paste('La altitud mínima es:',min(Altitud, na.rm = TRUE)))

- Crea una nueva df de datos sólo con las muestras de la raza Olotillo:
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
### Ejercicio 10: 
El resultado del Ejercicio 10 se encuentra en [Wiki sesion 1.3] (https://github.com/joaquintorresnunez/Tareas_BioinfRepro2026_JITN/wiki/Notas-sesion-1.3)


## Conclusion
En esta sesión se consolidaron las bases de la programación en R y el análisis genómico reproducible a través del desarrollo de scripts, la creación de funciones personalizadas y la manipulación de datos con dplyr. Asimismo, se abordaron metodologías clave en genética de poblaciones mediante la ejecución e interpretación de pruebas de Mantel y Mantel Parcial utilizando el paquete vegan. En conjunto, la resolución de estos ejercicios integró el procesamiento de datos biológicos con el uso de código anotado, una estructura clara de trabajo y la documentación en GitHub, elementos esenciales para garantizar la transparencia y la reproducibilidad científica. Dentro de los mayores desafíos al realizar dicha tarea destacaron la interpretación y comprensión de pruebas estadísticas asociadas al análisis de poblaciones, en conjunto con los ciclos de iteración; esto demandó un análisis más profundo del código utilizado, lo que a su vez permitió un mejor entendimiento de scripts más complejos.



