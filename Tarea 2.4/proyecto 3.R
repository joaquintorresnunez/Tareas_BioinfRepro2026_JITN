# ==============================================================================
# PROYECTO 3: Peopling of the Americas
# Analisis de olas migratorias y afinidades de poblacion mediante ADMIXTOOLS2
# ==============================================================================

library(admixtools)
install.packages("dplyr")
library(dplyr)
setwd("/home/joaquin/Escritorio/genantiguo")
library(admixtools)

# 1. Rutas a los genotipos crudos
geno_prefix <- "v62.0_1240k_public"
outdir_p3 <- "f2_project3"  # Carpeta nueva para este proyecto


# 2. Definir todas las poblaciones que necesitas
todas_las_pops <- c("Pima.DG", "Karitiana.DG", "Surui.DG", "PEL.DG", "MXL.DG", "USA_Anzick_realigned.SG", "USA_Ancient_Beringian.SG", "USA_Nevada_SpiritCave_11000BP.SG","Mbuti.DG", "CHB.DG", "Papuan.DG", "Denisova.DG")

# 3. Precomputar f2 para estas poblaciones específicas
extract_f2(
  pref      = geno_prefix,
  outdir    = outdir_p3,
  pops      = todas_las_pops,
  overwrite = TRUE
)

# 4. Cargar los bloques recién generados
f2_blocks <- f2_from_precomp(outdir_p3)

# 2. Definición de conjuntos de poblaciones
target3   <- c("Pima.DG", "Karitiana.DG", "Surui.DG", "PEL.DG", "MXL.DG")
source3   <- c("USA_Anzick_realigned.SG", "USA_Ancient_Beringian.SG", "USA_Nevada_SpiritCave_11000BP.SG")
outgroup3 <- c("Mbuti.DG", "CHB.DG", "Papuan.DG", "Denisova.DG")

# 3. Outgroup-f3: Deriva genética compartida entre targets y fuentes
# Medimos f3(Mbuti; Target, Source)
f3_results <- f3(
  data = f2_blocks,
  pop1 = "Mbuti.DG",
  pop2 = target3,
  pop3 = source3
)
print(f3_results)
# 4. Pruebas f4: Controles de asimetría y señal Australo-Melanesia
# 4a. Evaluacion de señal de afinidad con Oceanía (Papuan vs Pima)
f4_australasian <- f4(
  data = f2_blocks,
  pop1 = target3,
  pop2 = source3,
  pop3 = "Papuan.DG",
  pop4 = "Mbuti.DG"
)
print(f4_australasian)
# 4b. Evaluacion de afinidad diferencial entre fuentes antiguas respecto a targets
f4_sources <- f4(
  data = f2_blocks,
  pop1 = target3,
  pop2 = source3,
  pop3 = source3, 
  pop4 = "Mbuti.DG"  
)
print(f4_sources)
# 5. qpWave: Determinación de la cantidad de flujos de ancestría (rank)
wave_results <- list()
for (tgt in target3) {
  wave_results[[tgt]] <- qpwave(
    data = f2_blocks,
    left = c(tgt, source3),
    right = outgroup3
  )
}

for (tgt in target3) {
  cat("\n========================================\n")
  cat("Resultados qpWave para:", tgt, "\n")
  cat("========================================\n")
  print(wave_results[[tgt]])
}

# 6. qpAdm: Modelado de mezcla de 3 vías

adm_3way_results <- list()

for (tgt in target3) {
  adm_3way_results[[tgt]] <- qpadm(
    data   = f2_blocks,
    left   = source3,      # Fuentes antiguas
    right  = outgroup3,    # Outgroups
    target = tgt           # Target actual de la iteración
  )
}
for (tgt in names(adm_3way_results)) {
  cat("\n==================================================\n")
  cat(" MODELO qpAdm DE 3 VÍAS PARA:", tgt, "\n")
  cat("==================================================\n")
  
  cat("\n[1] Proporciones de mezcla (weights):\n")
  print(adm_3way_results[[tgt]]$weights)
  
  cat("\n[2] Ajuste del modelo y subtipos (popdrop):\n")
  print(adm_3way_results[[tgt]]$popdrop)
}
