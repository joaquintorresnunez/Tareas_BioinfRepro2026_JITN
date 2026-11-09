library(admixtools)
install.packages("dplyr")
library(dplyr)

setwd("/home/joaquin/Escritorio/genantiguo")

prefix <- "v62.0_1240k_public"
outdir_p2 <- "f2_project2"  # Carpeta nueva para este proyecto

# 1. Definir las poblaciones necesarias
target2   <- c("Russia_MLBA_Sintashta.AG", "Kazakhstan_Maitan_MLBA_Alakul.AG", "Russia_LBA_Srubnaya_Alakul.SG")
source2   <- c("Iran_GanjDareh_N.AG", "Russia_Sidelkino_HG.SG")
outgroup2 <- c("Mbuti.DG", "CHB.DG", "Papuan.DG", "Russia_UstIshim_IUP.DG", "Denisova.DG")

all_pops <- unique(c(target2, source2, outgroup2))

# 2. Extraer los bloques f2
extract_f2(
  pref = prefix,
  outdir = outdir_p2,
  pops = all_pops,
  transitions = FALSE,
  transversions = TRUE,
  overwrite = TRUE,
  blgsize = 0.05,
  verbose = TRUE
)

# 3. Cargar los bloques precalculados
f2_blocks <- f2_from_precomp(outdir_p2)

# 4. Calcular f3
f3_results <- f3(
  data = f2_blocks,
  pop1 = "Mbuti.DG",
  pop2 = target2,
  pop3 = source2
)
print(f3_results)

ref_pops <- c("Iran_GanjDareh_N.AG", 
              "Russia_Sidelkino_HG.SG")
# 5. Calcular f4
f4_results <- f4(
  data = f2_blocks,
  pop1 = target2,
  pop2 = ref_pops,
  pop3 = ref_pops,
  pop4 = "Mbuti.DG"
)
print(f4_results)

# 6. Evaluar rango de flujos con qpWave
wave2 <- qpwave(
  data = f2_blocks,
  left = c(target2[1], source2),
  right = outgroup2
)
print(wave2)

install.packages("dplyr")
library(dplyr)
# 7. Modelar la mezcla con qpAdm
qpadm_results <- map_dfr(target2, function(tg) {
  res <- qpadm(
    data = f2_blocks,
    left = source2,
    right = outgroup2,
    target = tg
  )
  
  res$weights %>% 
    mutate(
      target = tg,
      pvalue = res$popdrop$p[1]
    )
})

print(qpadm_results)
