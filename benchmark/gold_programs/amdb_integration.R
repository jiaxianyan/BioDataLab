# Install devtools if not available
# install.packages("remotes", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

# Install traitdata package from Github
# remotes::install_github("RS-eco/traitdata", build_vignettes = T, force=T)

library(traitdata)
data(mammal_diet2)

mus_trophic <- mammal_diet2$TrophicLevel[mammal_diet2$scientificNameStd == "Mus musculus"]
unique(mus_trophic)


mus_trophic <- mammal_diet2$TrophicLevel[mammal_diet2$scientificNameStd == "Propithecus verreauxi"]
unique(mus_trophic)
