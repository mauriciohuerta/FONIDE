rm(list = ls())

### Función para extraer los ruts desde los normbres de los archivos

nameling <- function(name){
  new.name <- gsub("Participant_", "P", name)
  new.name <- gsub("_2-S1-eRaw.RData", "_K", new.name)
  new.name <- gsub("_1-S1-eRaw.RData", "_M", new.name)
  new.name <- gsub("1-S1-eRaw.RData", "_M", new.name)
  new.name <- gsub("2-S1-eRaw.RData", "_K", new.name)
  return(new.name)
}

### Librerías a usar

library(eyelinker)
library(zoo)
library(dplR)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(grid)
library(data.table)
library(imputeTS)
#library(beepr)

### Establecer carpeta principal

setwd("C:/Users/mauri/Dropbox/Solicitudes de archivos/BDD/News")
main.path <- getwd()

### Carpetas secundarias

folders <- c("/C1/Kanjis", "/C1/Matematicas", "/C2/Kanjis", "/C2/Matematicas", "/C3/Kanjis", "/C3/Matematicas")

paths <- paste(main.path, folders, sep="")

### Cargando en memoria los archivos

for(f in paths){
  for(j in dir(f)){
    load(paste(f,j,sep="/"))
    assign(nameling(j), data.clean)
    rm(data.clean)
  }
}



sujetos <- ls()[!(ls() %in% c("f", "folders", "j", "main.path", "paths", "nameling"))]

ruts <- read.csv2("../00Backups_ArchivosRoot/correctas.csv")
ruts <- subset(ruts, select = c(archivo, Grupo))
ruts <- unique(ruts)
ruts$archivo <- gsub("Participant_", "", ruts$archivo)
ruts$archivo <- gsub("_2-S1-eRaw.csv", "", ruts$archivo)
ruts$archivo <- gsub("_2-S1-eRaw faltan 10 ult trials.csv", "", ruts$archivo)
ruts$archivo <- gsub("_1-S1-eRaw.csv", "", ruts$archivo)
ruts$archivo <- gsub("-S1-eRaw.csv", "", ruts$archivo)
ruts$archivo <- gsub("_S1-eRaw.csv", "", ruts$archivo)
ruts <- unique(ruts)
ruts$archivo <- as.character(ruts$archivo)

id9 <- which(nchar(ruts$archivo) == 9)
ruts$archivo[id9] <- gsub(".$", "", ruts$archivo[id9])
id10 <- which(nchar(ruts$archivo) == 10)
ruts$archivo[id10] <- gsub(".$", "", ruts$archivo[id10])
ruts$archivo[id10] <- gsub(".$", "", ruts$archivo[id10])
table(nchar(ruts$archivo))

for(s in sujetos){
  aux <- get(s)
  aux$subject_id <- as.character(unique(na.omit(aux$subject_id)))
  if(nchar(unique(aux$subject_id)) == 9){
    aux$subject_id <- gsub(".$", "", aux$subject_id)
  }
  if(nchar(unique(aux$subject_id)) == 10){
    aux$subject_id <- gsub(".$", "", aux$subject_id)
    aux$subject_id <- gsub(".$", "", aux$subject_id)
  }
  assign(s, aux)
  rm(aux)
}

names(ruts) <- c("subject_id", "Grupo")

for(s in sujetos){
  if(!(unique(get(s)$subject_id) %in% ruts$subject_id)){
    cat("CHEK ",s)
  }
  aux <- get(s)
  aux <- merge(aux, ruts, by="subject_id", all.x=T, all.y=F)
  assign(s, aux)
}

rm(ruts, id9, id10)

# Ajustando el sampling

relleno <- data.frame(relative_time = 1:1050, var=NA)

for(s in sujetos){
  aux <- subset(get(s), select = c(subject_id, tipo_prueba,relative_time,Pupil.Diameter.Left,Pupil.Diameter.Right,trial_nr, Grupo))
  laux <- list()
  for(j in 1:60){
    cat(s, j, "\n")
    aux.i <- subset(aux, trial_nr == j)
    aux.i$relative_time <- floor(aux.i$relative_time/10)
    index <- !duplicated(aux.i$relative_time)
    aux.i <- aux.i[index,]
    DPaste <- merge(aux.i, relleno, by = "relative_time", all.x=T, all.y=T)
    DPaste$subject_id <- unique(na.omit(DPaste$subject_id))
    DPaste$tipo_prueba <- unique(na.omit(DPaste$tipo_prueba))
    DPaste$trial_nr <- unique(na.omit(DPaste$trial_nr))
    end <- which(!is.na(DPaste$Pupil.Diameter.Left))
    if(length(end) > 0){
      DPaste$Pupil.Diameter.Left[1] <- DPaste$Pupil.Diameter.Left[end[1]]
      end <- end[length(end)]
      DPaste <- DPaste[1:end,]
      DPaste$ps <- na_interpolation(DPaste$Pupil.Diameter.Left, option = "spline")
    } else{
      DPaste$ps <- DPaste$Pupil.Diameter.Left
    }
    laux[[j]] <- DPaste
    rm(aux.i, index, DPaste, end)
  }
  assign(paste(s,"sample",sep="."), subset(rbindlist(laux), select = -var))
  rm(aux, laux)
  rm(list = ls()[ls() %in% s])
}

rm(relleno)

sujetos <- paste(sujetos, "sample",sep=".")

### Suavizamiento y normalizacion Z

for(s in sujetos){
  aux <- get(s)
  laux <- list()
  for(j in 1:60){
   # cat(s, j, "\n")
    aux.i <- subset(aux, trial_nr == j)
    if(sum(is.na(aux.i$ps)) == nrow(aux.i)){
      cat("CHEK ",s,j,"\n")
      aux.i$FRps <- NA
      aux.i$Zps <- NA
    } else{
    aux.i$FRps <- detrend.series(aux.i$ps, method="Friedman", return.info=T, make.plot = F)$curves
    aux.i$Zps <- scale(aux.i$FRps, center = TRUE, scale = TRUE)
    }
    laux[[j]] <- aux.i
    rm(aux.i)
  }
  assign(paste(s,"Z",sep="_"), rbindlist(laux))
  rm(aux, laux)
  rm(list = ls()[ls() %in% s])
}

rm(P209147572_K.sample_Z,  P212693634_K.sample_Z, P215795683_K.sample_Z, P21802476_M.sample_Z, P21820393_M.sample_Z,  P258568575_K.sample_Z)

sujetos <- paste(sujetos,"Z",sep="_")
stdrop <- c("P209147572_K.sample_Z",  "P212693634_K.sample_Z", "P215795683_K.sample_Z", "P21802476_M.sample_Z", "P21820393_M.sample_Z",  "P258568575_K.sample_Z")

id <- which(sujetos %in% stdrop)
sujetos <- sujetos[-id]

rm(id, stdrop)
# uniendo todas las bases de datos

laux <- list()
for(j in 1:length(sujetos)){
  print(j)
  laux[[j]] <- get(sujetos[j])
}
full_data <- rbindlist(laux)
rm(laux)

full_kanjis <- subset(full_data, tipo_prueba == "kanjis")
full_maths <- subset(full_data, tipo_prueba == "matematicas")

kanjis.C <- subset(full_kanjis, Grupo == "control")
kanjis.T <- subset(full_kanjis, Grupo == "treatment")
math.C <- subset(full_maths, Grupo == "control")
math.T <- subset(full_maths, Grupo == "treatment")

# Funcion que ajusta a un baseline especifico todos los datos
  
baseline <- function(baseline, data){
  subjects <- unique(data$subject_id)
  lauxGen <- list()
  for(i in 1:length(subjects)){
    aux <- subset(data, subject_id == subjects[i])
    laux <- list()
    for(j in 1:60){
    #  cat("sujeto ",i, " de ", length(subjects), " trial ", j, "\n")
      aux.i <- subset(aux, trial_nr == j)
      if(sum(is.na(aux.i$ps)) == nrow(aux.i)){
        aux.i$ZpsBL <- NA
      } else{
        aux.i$ZpsBL <- scale(aux.i$Zps, center = aux.i$Zps[baseline], scale = TRUE)
      }
      laux[[j]] <- aux.i
      rm(aux.i)
    }
    lauxGen[[i]] <- rbindlist(laux)
    rm(aux, laux)
  }
  BLdata <- rbindlist(lauxGen)
  return(BLdata)
}

# Funcion que calcula promedios
promedios <- function(datos){
  idlimit <- which(datos$relative_time < 1050)
  datos <- datos[idlimit,]
  subjects <- unique(datos$subject_id)
  laux <- list()
  for(i in 1:length(subjects)){
    aux <- subset(datos, subject_id == subjects[i])
    medias.i <- aggregate(aux$ZpsBL, by=list(aux$relative_time),
                                FUN = "mean", na.action = NULL, na.rm=T)
    laux[[i]] <- medias.i
    rm(aux, medias.i)
  }
  meansSujetos <- rbindlist(laux)
  names(meansSujetos) <- c("relative_time", "ZpsBL")
  rm(laux)
  MediaGeneral <- try(aggregate(meansSujetos$ZpsBL, by=list(meansSujetos$relative_time),
                                FUN = "mean", na.action = NULL, na.rm=T))
  SDGeneral <- try(aggregate(meansSujetos$ZpsBL, by=list(meansSujetos$relative_time),
                                FUN = "sd", na.rm=T))
  names(MediaGeneral) <- c("relative_time", "media")
  names(SDGeneral) <- c("relative_time", "SD")
  outline <- merge(MediaGeneral, SDGeneral, by="relative_time", all=T)
  return(outline)
}
  

# Función que calcula diferencias

maxdif <- function(Control, Tratamiento, BaseLines){
  aux <- list()
  i <- 0
  for(bl in BaseLines){
    cat(bl, "\n")
    i <- i + 1
    BLC <- baseline(bl, data = Control)
    BLT <- baseline(bl, data = Tratamiento)
    Co <- promedios(BLC)
    Tr <- promedios(BLT); print(Tr$SD)
    pC <- Co$media
    pT <- Tr$media    
    isC <- pC +Co$SD
    iiC <- pC -Co$SD
    isT <- pT +Tr$SD
    iiT <- pT -Tr$SD
    len <- min(length(pC), length(pT))
    idMaxDif <- which.max(abs(pC[1:len] - pT[1:len]))
    aux[[i]] <- data.frame(BL = bl, PC = pC[idMaxDif], PT = pT[idMaxDif], IIC = iiC[idMaxDif], ISC = isC[idMaxDif],
                           IIT = iiT[idMaxDif], IST = isT[idMaxDif], where = idMaxDif)
  }
  output <- rbindlist(aux)
  return(output)
}

maxdif.int <- function(Control, Tratamiento, BaseLines){
  aux <- list()
  i <- 0
  for(bl in BaseLines){
    cat(bl, "\n")
    i <- i + 1
    BLC <- baseline(bl, data = Control)
    BLT <- baseline(bl, data = Tratamiento)
    Co <- promedios(BLC)
    Tr <- promedios(BLT); print(Tr$SD)
    pC <- Co$media
    pT <- Tr$media    
    isC <- pC + Co$SD
    iiC <- pC - Co$SD
    isT <- pT + Tr$SD
    iiT <- pT - Tr$SD
    len <- min(length(pC), length(pT))
    idMaxDif <- which.max(abs(pC[1:len] - pT[1:len]))
    aux[[i]] <- data.frame(BL = bl, PC = pC[idMaxDif], PT = pT[idMaxDif], IIC = iiC[idMaxDif], ISC = isC[idMaxDif],
                           IIT = iiT[idMaxDif], IST = isT[idMaxDif], where = idMaxDif)
  }
  output <- rbindlist(aux)
  return(output)
}

sal <- maxdif(Control = kanjis.C, Tratamiento = kanjis.T, BaseLines = c(1, seq(10, 1040, by = 10)))

sal

sal$cte <- (sal$PC + sal$PT)/2

plot(sal$BL, sal$PC - sal$cte, type="l", ylim=c(-3, 3))
lines(sal$BL, sal$PT - sal$cte)
abline(h=0, col=2)

ggplot() + #ylim(-2,2) +
  geom_line(aes(x=sal$BL/100, y=sal$PC - sal$cte), color = "red") +
  geom_line(aes(x=sal$BL/100, y=sal$IIC/sqrt(104) - sal$cte), color = "red", alpha=0.2) + 
  geom_line(aes(x=sal$BL/100, y=sal$ISC/sqrt(104) - sal$cte), color = "red", alpha=0.2) +
  geom_line(aes(x=sal$BL/100, y=sal$PT - sal$cte), color = "blue") +
  geom_line(aes(x=sal$BL/100, y=sal$IIT/sqrt(116) - sal$cte), color = "blue", alpha=0.2) + 
  geom_line(aes(x=sal$BL/100, y=sal$IST/sqrt(116) - sal$cte), color = "blue", alpha=0.2) + 
  geom_ribbon(aes(x= sal$BL/100, ymin = sal$IIC/sqrt(104) - sal$cte, ymax = sal$ISC/sqrt(104) - sal$cte), fill = "red", alpha=0.2) +
  geom_ribbon(aes(x= sal$BL/100, ymin = sal$IIT/sqrt(116) - sal$cte, ymax = sal$IST/sqrt(116) - sal$cte), fill = "blue", alpha=0.2) +
  geom_hline(yintercept = 0)

ggplot() + #ylim(-2,2) +
  geom_line(aes(x=sal$BL/100, y=sal$PC), color = "red") +
  geom_line(aes(x=sal$BL/100, y=sal$IIC), color = "red", alpha=0.2) + 
  geom_line(aes(x=sal$BL/100, y=sal$ISC), color = "red", alpha=0.2) +
  geom_line(aes(x=sal$BL/100, y=sal$PT), color = "blue") +
  geom_line(aes(x=sal$BL/100, y=sal$IIT), color = "blue", alpha=0.2) + 
  geom_line(aes(x=sal$BL/100, y=sal$IST), color = "blue", alpha=0.2) + 
  geom_ribbon(aes(x= sal$BL/100, ymin = sal$IIC, ymax = sal$ISC), fill = "red", alpha=0.2) +
  geom_ribbon(aes(x= sal$BL/100, ymin = sal$IIT, ymax = sal$IST), fill = "blue", alpha=0.2) +
  geom_hline(yintercept = 0)
  

  
ggplot() + 
  geom_line(aes(x=sal$BL/100, y=sal$PC), color="blue") +
  geom_point(aes(x=sal$BL/100, y=sal$PC), color="blue")+
  geom_errorbar(aes(x= sal$BL/100, ymin=sal$IIC, ymax=sal$ISC), width=.1, col="blue", alpha = 32) +
  geom_line(aes(x=sal$BL/100, y=sal$PT), color="red") +
  geom_point(aes(x=sal$BL/100, y=sal$PT), color="red")+
  geom_errorbar(aes(x= sal$BL/100, ymin=sal$IIT, ymax=sal$IST), width=.1, col="red", alpha=.3)+
  xlab("Baseline (time in seconds)") + 
  ylab("Maximum difference (Z-Score)")
