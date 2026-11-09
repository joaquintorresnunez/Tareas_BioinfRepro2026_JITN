###Cargar librerias, problemas con cargar tidyverse, se decide solo cargar dplr##
setwd("/home/joaquin/Escritorio/genantiguo/aadr_1000G_f2_all")
library(admixtools)
install.packages("dplyr")
library(dplyr)

####Cargar f2#####
outdir <- "/home/joaquin/Escritorio/genantiguo/aadr_1000G_f2_all"
f2_blocks <- f2_from_precomp(outdir)

####Definir poblaciones de referencia####
source1 <- c("Turkey_Marmara_Barcin_N.AG", 
             "Russia_Samara_EBA_Yamnaya.AG", 
             "Luxembourg_Mesolithic.DG")

####Definir poblaciones Outgroups#### 
outgroup1 <- c("Mbuti.DG", "CHB.DG", "Papuan.DG", "Russia_UstIshim_IUP.DG", "Denisova.DG")

####Definir poblaciones target####
targets <- c("GBR.DG", "IBS.DG", "TSI.DG", "FIN.DG", "Sardinian.DG")

####Estimacion de la deriva genica compartida (estadistico f3)####
# Calculamos f3 entre GBR.DG y las 3 fuentes utilizando Mbuti.DG como Outgroup
f3_results <- f3(f2_blocks, pop1 = "Mbuti.DG", pop2 = "GBR.DG", pop3 = source1)

print("Resultados de f3 para GBR.DG:")
print(f3_results)

####Definición de referencias ancestrales para iterar e calculo f4####
ref_pops <- c("Luxembourg_Mesolithic.DG", 
              "Turkey_Marmara_Barcin_N.AG", 
              "Russia_Samara_EBA_Yamnaya.AG")

####Genera todas las combinaciones de poblaciones de referencia f4 (Mbuti, GBR; Ref_A, Ref_B)####
f4_results <- f4(f2_blocks,
                  pop1 = "GBR.DG",
                  pop2 = ref_pops,
                  pop3 = ref_pops,
                  pop4 = "Mbuti.DG")

print("Resultados de f4:")
print(f4_results)

####Evaluación de rango con qpwave####
wave1 <- qpwave(f2_blocks,
                left = c("GBR.DG", source1),
                right = outgroup1)

print("Resultados de qpWave:")
print(wave1)

####Comparación de modelos qpAdm (2-way vs 3-way en GBR.DG)####

# Modelo de 2 fuentes (Farmer + Steppe)
admix_2way <- qpadm(f2_blocks, 
                    left = source1[1:2], 
                    right = outgroup1, 
                    target = "GBR.DG")

# Modelo de 3 fuentes (Farmer + Steppe + WHG)
admix_3way <- qpadm(f2_blocks, 
                    left = source1, 
                    right = outgroup1, 
                    target = "GBR.DG")

print(admix_2way$weights)
print(admix_3way$weights)

#### Iteración enfocado únicamente en la comparación de Admixture (2-way vs 3-way) ####
for (target_pop in targets) {
  
  cat("\n=======================================================\n")
  cat(" COMPARACIÓN DE MEZCLA (ADMIXTURE) PARA:", target_pop, "\n")
  cat("=======================================================\n")
  
  # Cálculo de modelo 2-way
  adm_2 <- qpadm(f2_blocks, 
                 left = source1[1:2], 
                 right = outgroup1, 
                 target = target_pop)
  
  # Cálculo de modelo 3-way
  adm_3 <- qpadm(f2_blocks, 
                 left = source1, 
                 right = outgroup1, 
                 target = target_pop)
  
  # Impresión exclusiva de las tablas de pesos ($weights)
  cat("\n--- Pesos Modelo 2-way ---\n")
  print(adm_2$weights)
  
  cat("\n--- Pesos Modelo 3-way ---\n")
  print(adm_3$weights)
  cat("\n")
}


https://www.youtube.com/watch?v=jjKJsOtZ5Vk&t=934s

