##préalables avant exécution du script
#installer R + Rstudio ici https://rstudio.com/products/rstudio/download/
#lancer et installer le package Fragman avec install.packages("Fragman")

##script
#Thibault CARON, 17/06/2020
#lire et analyser les fichier de génotypage .fsa en alternative à GeneMapper (ici un exemple pour 11 microsatellites en quadruplex)

library(Fragman) #charger le package
help("Fragman") #documentation du package
setwd("/media/thibault/TOSHIBA/supplements_croq/QTL/INDEL/1.raw data/E149_Envoi_du_20180314/E149_Envoi_du_20180314/E149-Plaque1-20180313/C8") #changer le répertoire de travail (où se trouvent les .fsa)
#travailler de préférence sur disque dur interne (plus rapide)

#1. lire les .fsa : folder doit être le même répertoire que setwd(), soit "."
fsa <- storing.inds(folder = ".", channels = 5, fourier = T, saturated = T, lets.pullup = T, plotting = T, rawPlot = F, llength = 3000,
                    ulength = 80000)

#2. renseigner l'échelle
ladder <- c(50,75,100,139,150,160,200,250,300,340,350,400,450,490,500) #tailles des pics
ladder.info.attach(fsa, ladder = ladder, channel.ladder = 5, method = "iter2", ladd.init.thresh = NULL, env = parent.frame(), prog = T, draw = T, 
                   attempt = 5) #attacher l'échelle à chaque fichier .fsa, ne pas tenir compte de l'éventuel message d'avertissement

#3. visualiser les génotypes des parents ou d'un ensemble d'individus pour un fluorophore spécifique pour préparer le choix des panels (sélectionner une seule chaîne d'un seul multiplex à la fois)
overview2(my.inds = fsa[substr(names(fsa),2,2)==1], channel = 1, ladder = ladder, channel.ladder = 5, xlim = c(100,250), ploidy = 1, method = "cor")
overview2(my.inds = fsa[23], channel = 3, ladder = ladder, channel.ladder = 5, xlim = c(100,250), ploidy = 1, method = "cor", lwd = 2)

#3.5 permet d'enregistrer la taille du pic attendu dans l'objet de son choix, après apparition de la croix, cliquer puis echap ou finish
#info$size?[?] <- locator(type="p", pch=20, col="red")$x 

#4. renseigner les tailles et poids attendus ? (renseigner ici)
#info <- read.csv("info.csv",sep="\t") #(possibilité d'importer un tableur .csv dans le répertoire de travail avec les mêmes noms de champs avec read.csv())
info <- data.frame("multiplex" = c(rep(1,4),rep(2,4),rep(3,4)),
                        "marker" = c("Proq12","Proq13","Proq80","Proq83","Proq73","Proq75","Proq77","Proq81","Proq74","Proq78","Proq88","Proq93"),
                        "dye" = c("G","Y","R","B","B","Y","G","R","Y","G","B","R"),
                        "pattern" = c("AT","AT","CTG","ATG","TTC","TTC","CTC","CTG","TTC","ACT","AAC","AAG"),
                        "size1" = c(183, 185, 175, 175, 179, 143, 195, 177, 157, 199, 155, 165),
                        "size2" = c(185, 189, 175, 181, 187, 162, 195, 183, 163, 199, 182, 177))
info$marker <- as.vector(info$marker) #colonne des marqueurs en vecteur (pour manipuler ensuite)
channels <- c("B"=1, "G"=2, "Y"=3, "R"=4) #chaque chaine correspond à un fluorophore : 1 = blue, 2 = green, 3 = yellow, 4 = red, 5 = orange (ladder)
info$channel <- channels[as.vector(info$dye)] #ajouter les numéros de chaînes correspondant aux fluorophores

#4. mesure automatisée des tailles pour chaque chaîne de chaque multiplex
results <- data.frame() #initialiser le tableau des résultats
log <- NULL #initialiser le fichier de rapport
for(multiplex in 1:3){ #pour chaque multiplex de 1 à 3
#multiplex <- 1 #débugage boucle
  for (channel in 1:4){ #pour chaque chaîne de 1 à 4
  #channel <- 1 #débugage boucle
    name <- info$marker[info$multiplex==multiplex & info$channel==channel] #récupérer le nom du marqueur
    log <- c(log, print(paste("traitement du marqueur ", name, " en multiplex ", multiplex," chaîne ", channel, " / ", info$dye[info$multiplex==multiplex & info$channel==channel], " : ", sep = "")))
    #récupérer les noms de souches, à adapter selon le format des noms de fichiers (ici séparés par "_", noms en 3ème position)
    strains <- strsplit(names(fsa)[1:length(names(fsa))/3], "_") #séparation du champ de noms par "_" en listes
    strains <- unique(unlist(lapply(strains,function(x){x[[3]]}))) #récupérer le troisième élément de chaque liste
    short <- info$size1[info$channel==channel & info$multiplex==multiplex] #soit "short" la petite taille attendue
    long <- info$size2[info$channel==channel & info$multiplex==multiplex] #soit "long" la grande taille attendue
    log <- c(log, print(paste("petite et grande taille attendues : ", short, " et ", long, sep = "")))
    #définir un premier panel large à partir des tailles attendues
    #on construit une séquence à partir de 5 base en dessous de la petite taille attendue, jusqu'à 5 bases au dessus de la grande taille attendue,
    #avec un pas de 0.1 base
    panel1 <- seq(short-5, long+5, 0.01)
    log <- c(log, print(paste("construction du premier panel de ", short-5, " à ", long+5, " tous les 0.01 bases", sep = "")))
    #recherche des scores : premier passage large
    #1er argument + n.inds : les noms des fichiers doivent correspondre à la sélection du multiplex (ici notés "M1etc")
    #en ploidy 1 ne recherche qu'un seul pic (inutile de paramétrer left.cond, right.cond et shift)
    #window : taux d'erreur en base // panel, laisser une faible valeur
    #plotting : T/F pour afficher les électrogrammes traités (!= electro : électrogrammes bruts)
    scores1 <- score.markers(fsa[substr(names(fsa),2,2)==multiplex], channel = channel, n.inds = c(1:length(fsa[substr(names(fsa),2,2)==multiplex])),
                            panel = panel1, ladder = ladder, channel.ladder = 5, ploidy = 1, shift = 1.5,  left.cond = c(0.3,5), right.cond = 0.5,
                            warn = F, window = 0.1, init.thresh = 5000, ladd.init.thresh = NULL, method = "iter2", env = parent.frame(), 
                            my.palette = NULL, plotting = T,  electro = F, pref = 3)
    #contruire de nouveaux panels à partir des tailles observées
    observed <- unlist(lapply(scores1, function(x){x$wei})) #extraire les tailles observées
    #distrib <- as.data.frame(table(as.integer(observed))) #calculer la distribution
    distrib <- as.data.frame(table((round(observed/0.5))*0.5)) #calculer la distribution en arrondissant à 0.5 base
    colnames(distrib)[1]<-"taille"
    log <- c(log, print("premier passage effectué, la distribution des pics est la suivante : taille / fréquence"))
    dis <- apply(distrib,1,function(x)paste(x[1],"\t",x[2],sep=""))
    log <- c(log, print(paste(dis, sep = "")))
    if(dim(distrib)[1]>=2){ #s'il y a au moins 2 tailles de pic obtenues
      #ordonner par fréquences décroissantes et sélectionner les 2 plus fréquents
      distrib2 <- distrib[order(distrib[,2],decreasing = T),][1:2,]
      log <- c(log, print(paste("les deux tailles plus fréquemment observées sont de : ", distrib2[1,1], " et ", distrib2[2,1], " bases", sep = "")))
      #contruire un deuxième panel plus précis en fonction des tailles observées
      #on prend de 8 bases en dessous et jusqu'à 1 base au dessus du petit pic, puis de 1 base en dessous et jusqu'à 5 bases au dessus du grand pic, avec un pas de 1 base
      panel2 <- c(seq(min(as.numeric(as.vector(distrib2[,1])))-2, min(as.numeric(as.vector(distrib2[,1])))+2,1),
                  seq(max(as.numeric(as.vector(distrib2[,1])))-2, max(as.numeric(as.vector(distrib2[,1])))+2,1))
      log<-c (log, print(paste("construction du panel définitif de ", min(as.numeric(as.vector(distrib2[,1])))-2, " à ", min(as.numeric(as.vector(distrib2[,1])))+2,
                  " et de ", max(as.numeric(as.vector(distrib2[,1])))-2, " à ", max(as.numeric(as.vector(distrib2[,1])))+2, " tous les 1 bases", sep = "")))}
    if(dim(distrib)[1]==1){ #s'il n'y a qu'une seule taille de pic obtenue
      log <-c (log, print(paste("l'unique taille observée est de : ", distrib[1,1], " bases", sep = "")))
      #on prend de 5 bases en desous jusqu'à 5 bases au dessus du pic unique, avec un pas de 1
      panel2 <- seq(as.numeric(as.vector(distrib[1,1]))-5, as.numeric(as.vector(distrib[1,1]))+5,1)
      log <- c(log, print(paste("construction du panel définitif de ", as.numeric(as.vector(distrib[1,1]))-5, " à ", as.numeric(as.vector(distrib[1,1]))+5, " tous les 1 bases", sep = "")))}
    #deuxième recherche des scores plus précise : window à 1 base + ploidy = 2 (certains cas où 2 pics sont clairement visibles aux tailles observées)
    scores2 <- score.markers(fsa[substr(names(fsa),2,2)==multiplex], channel = channel, n.inds = c(1:length(fsa[substr(names(fsa),2,2)==multiplex])),
                             panel = panel2, ladder = ladder, channel.ladder = 5, ploidy = 2, shift = 1.5,  left.cond = c(0.3,5), right.cond = 0.3,
                             warn = F, window = 1, init.thresh = 5000, ladd.init.thresh = NULL, method = "iter2", env = parent.frame(), 
                             my.palette = NULL, plotting = T,  electro = F, pref = 3)
    geno <- get.scores(scores2) #récupérer le tableau de génotypes
    log <- c(log, print("deuxième passage effectué, élimination des éventuels cas de 2 pics différents et remplacement des tailles obtenues par les tailles attendues"))
    #fusionner les 2 pics récupérés en 1 (sauf si tailles différentes -> NA)
    geno <- as.data.frame(apply(geno,1,function(x){ifelse(x[1]==x[2],x[1],NA)}))
    #remplacer les tailles obtenues par les tailles attendues
    geno2 <- NULL
    for(i in 1:dim(geno)[1]){ #pour tous les pics du marqueur
      if(!is.na(geno[i,1])){ #si le génotype n'est pas NA
        if(geno[i,1]<= short){ #si la taille observée est inférieure ou égale à la petite taille attendue
          geno2<-c(geno2,short) #remplacer par la petite taille attendue
          log<-c(log,print(paste("individu ", strains[i], " (n°", i, ") : taille obtenue de ", geno[i,1], " bases : remplacement par la petite taille attendue de ", short, " bases", sep = "")))}
        if(geno[i,1]>= long){ #si la taille observée est supérieure ou égale à la grande taille attendue
          geno2<-c(geno2,long) #remplacer par la grande taille attendue
          log <- c(log, print(paste("individu ", strains[i], " (n°", i, ") : taille obtenue de ", geno[i,1], " bases : remplacement par la grande taille attendue de ", long, " bases", sep = "")))}
        if(geno[i,1]> short & geno[i,1] < long){ #si la taille obtenue est entre les deux valeurs attendues :
          diffG <- geno[i,1]-short #soit diffG l'écart à gauche
          diffD <- long-geno[i,1] #soit diffD l'écart à droite
          #ici on règle le seuil de rapprochemment du pic observé par rapport aux tailles attendues, soit la distance entre le pic obtenu et l'une des valeurs attendues pour confirmer son remplacement
          seuil<-2 #plus cette valeur est grande, plus on demande au pic obtenu d'être proche d'un des 2 pics attendus pour être remplacé, sinon on met NA (voir ci-dessous)
          if(diffG>=seuil*diffD){ #si l'écart à gauche est au moins "seuil" fois plus important que celui de droite
            geno2 <- c(geno2, long) #on prend donc la grande taille attendue
            log <- c(log, print(paste("individu ", strains[i], " (n°", i, ") : taille obtenue de ", geno[i,1], " bases : remplacement par la grande taille attendue de ", long, " bases", sep = "")))}
          if(diffD>=seuil*diffG){ #si l'écart à droite est au moins "seuil" fois plus important que celui de gauche
            geno2 <- c(geno2, short) #on prend donc la petite taille attendue
            log <- c(log, print(paste("individu ", strains[i], " (n°", i, ") : taille obtenue de ", geno[i,1], " bases : remplacement par la petite taille attendue de ", short, " bases", sep = "")))}
          if(diffG<seuil*diffD & diffD<seuil*diffG){ #si l'écart n'est pas "seuil" fois plus important à droite qu'à gauche
            geno2 <- c(geno2, NA) #trop d'incertitude, on remplace par NA
            log <- c(log, print(paste("individu ", strains[i], " (n°", i, ") : taille obtenue de ", geno[i,1], " bases : trop d'incertitude, remplacement par NA", sep = "")))} 
        }
      }
      else{ #si la taille est NA
        geno2<-c(geno2,NA) #on laisse en NA
        log<-c(log,print(paste("individu ", strains[i], " (n°", i, ") : pas de pic détecté (NA)", sep = "")))}
    }
    geno2 <- as.data.frame(geno2)
    colnames(geno2) <- name #renomer le nom du marqueur dans le tableau de génotypes
    rownames(geno2) <- strains #renomer le tableau de génotypes avec les noms de souches
    #remplir le tableau de résultats
    if(nrow(results)==0){ #si le nombre de lignes du tableau des résultats est égal à 0 (tableau vide, cas du premier marqueur)
      results <- geno2 #simplement copier la première colonne dans le tableau
    }
    else { #sinon (pour les autres marqueurs)
      results <- merge(results, geno2, by = "row.names", all.x = T) #fusionner les nouveaux génotypes au tableau des résultats par noms de souches
      rownames(results) <- results$Row.names #copier les noms de souches en noms de lignes
      results <- results[,-1] #(enlever la colonne inutile)
    }
  }
  if(multiplex==3 & channel==4){ #après le dernier passage
    write(log, file = "rapport.txt") #imprimer le rapport
    print("terminé ; le fichier 'rapport.txt' contient une copie du rapport de détection (ou entrez 'log')")}
}


#5. vérifier la bonne affiliation des pics aux tailles attendues
#permet de vérifier ou changer la taille du pic observé dans l'objet de son choix, après apparition de la croix, cliquer puis echap ou finish
locator(type="p", pch=20, col="red")$x 
    
print(paste("missing values", sum(apply(results,2,function(x){sum(is.na(x))}))/(dim(results)[1]*dim(results)[2]), sep = " : ")) #calculer le % de NA
write.table(results, file = "results.txt", sep = "\t") #imprimer les résultats (tabulé)
