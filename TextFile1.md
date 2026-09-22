# Tarea sesión 3

## Ejercicio1
````
log
x=log10(50)
y=5
z=x + y"./BioinfinvRepro-master/BioinfinvRepro-master/Unidad1/Sesion3/PracUni1Ses3/maices/meta/submatJT.csv")
print(z)
````
## Ejercicio 2
```
x=(1:150)
y=x+2
print(y)
```

## Ejercicio 3
```
?sum
vector=c(-13432:234)
sum(vector>20)
```

## Ejercicio 4
```
getwd()
?read.delim
# Cargar el archivo de metadatos
meta_maiz <- read.delim(""./BioinfinvRepro-master/BioinfinvRepro-master/Unidad1/Sesion3/PracUni1Ses3/maices/meta/maizteocintle_SNP50k_meta_extended.txt")
````

## Ejercico 5
### ejemplo
```
for (i in 2:10){
  print(paste(i, "elefantes se columpiaban sobre la tela de una araña"))
}
```
### Ejercicio 5.1
```
for (i in 1:10) {
  x = 35/i
  print(paste("El resultado de la división 35 /", i, "es igual a", x))
}
```

### Ejercicio 5.2
```
for (i in 1:10) {
  if (i %% 2 == 0) {
    next
  }
  x <- 35/i
  print(paste("El resultado de la división 35 /", i, "es igual a", x))
}
```

### Ejercicio 5.3
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
### Ejercicio 6
#### ¿qué hacen los dos for loops del script?: 
#### ¿qué paquetes necesitas para correr el script?: (ade4), (ggplop), (sp)
#### ¿qué archivos necesitas para correr el script?: (read.fst_summary_fix.R), (read.effdist.R), (surveyed_mountains.tsv), (BerSS.sumstats.tsv), (Balpina_focalpoints.txt), 

### Ejercicio 7
```
#Red Bananas
calc.tetha=function(Ne, u){
  tetha=4*Ne*u
  return(tetha)
}
```
### Ejercicio 8

#Elefante blanco

