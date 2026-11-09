library(admixtools)
library(dplyr)
setwd("/home/joaquin/Escritorio/genantiguo")

prefix <- "v62.0_1240k_public"
outdir_all <-"f2_project4.1"
# 1. Definir los conjuntos del Proyecto 4 (ajustando los nombres según la salida del Step 1)
target4    <- "IBS.DG"

sources4   <- c(
  "Morocco_Iberomaurusian.AG", 
  "YRI.DG"
)

outgroups4 <- c("Mbuti.DG", "Papuan.DG") 
# 2. Unificar todas las poblaciones únicas
todas_las_pops_p4 <- unique(c(target4, sources4, outgroups4))

# 3. Extraer bloques f2 sobrescribiendo el directorio previo
extract_f2(
  pref      = "v62.0_1240k_public",
  outdir    = "f2_project4.1",
  pops      = todas_las_pops_p4,
  maxmiss   = 1,
  auto_only = FALSE,
  overwrite = TRUE
)

# 4. Cargar el objeto f2 actualizado en memoria
f2_blocks_p4 <- f2_from_precomp("f2_project4.1")

# 5. Outgroup-f3: Deriva genética compartida entre targets y fuentes
# Medimos f3(Mbuti; Target, Source)
f3_results_p4 <- f3(
  data = f2_blocks_p4,
  pop1 = "Mbuti.DG",
  pop2 = target4,
  pop3 = sources4
)
print(f3_results_p4)


# 6. Calculo f4. Evaluacion de afinidad diferencial entre fuentes antiguas respecto a targets
f4_sources <- f4(
  data = f2_blocks_p4,
  pop1 = target4,
  pop2 = sources4,
  pop3 = sources4, 
  pop4 = "Mbuti.DG"  
)
print(f4_sources)

####7. Calculo qpwave####
wave_results_4 <- list()
for (tgt in target4) {
  wave_results_4[[tgt]] <- qpwave(
    data = f2_blocks_p4,
    left = c(tgt, sources4),
    right = outgroups4
  )
}

for (tgt in target4) {
  cat("\n========================================\n")
  cat("Resultados qpWave para:", tgt, "\n")
  cat("========================================\n")
  print(wave_results_4[[tgt]])
}


####8. Calculo qpadmin####

adm_results_4 <- list()

for (tgt in target4) {
  adm_results_4[[tgt]] <- qpadm(
    data   = f2_blocks_p4,
    left   = sources4,     # Fuentes de ancestralidad (Morocco_Iberomaurusian.AG, YRI.DG)
    right  = outgroups4,   # Panel de outgroups (Mbuti.DG, Papuan.DG)
    target = tgt           # Población objetivo (IBS.DG)
  )
}

# Visualizar resultados
for (tgt in names(adm_results_4)) {
  cat("\n==================================================\n")
  cat(" MODELO qpAdm PARA TARGET:", tgt, "\n")
  cat("==================================================\n")
  
  cat("\n[1] Proporciones de mezcla (weights):\n")
  print(adm_results_4[[tgt]]$weights)
  
  cat("\n[2] Evaluación de ajuste global y submodelos (popdrop):\n")
  print(adm_results_4[[tgt]]$popdrop)
}
