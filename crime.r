rm(list=ls()) # Efface les variables en m????moire
#########################################################

Crimes <- read.table("crimestatfinal.csv",header=TRUE, sep=",")
dim (Crimes)
# on realise une regression sur toutes nos variables explicatives
regression <- lm (homicide ~ death_penalty + white + black + female + pauvrete + revenu + chomage + alloc, Crimes)
#print(regression)
print(summary(regression))

cor(Crimes$alloc, Crimes$chomage)
cor(Crimes$alloc, Crimes$revenu)
cor(Crimes$alloc, Crimes$pauvrete)
# on supprime les variables non determinantes au seuil de 1% 
regression2 <- lm (homicide ~ revenu+ death_penalty + white + black + female + pauvrete + chomage, Crimes)
#print(regression2)
print(summary(regression2))


# test de kolmogorov par rapport a une N(26, 638) et par rapport a une chi2.
hist(Crimes$homicide,  main = 'Taux d homicide pour 10 000 hab par county')

normale <-rnorm(3144,mean=26,sd=25)
chi2 <- 100*rchisq(3144,1)
cauch <- rcauchy(3144, 0, 0.01)
ex=rexp(3144, 2)
histExp =hist(ex)

ks.test(Crimes[2], ex)
ks.test(Crimes[2], normale)
ks.test(Crimes[2], chi2)
hist(chi2)

# On veut savoir si le taux homicide suit la meme loi que le county ait la peine de mort
# ou non
OuiPeine =  Crimes$homicide[Crimes$death_penalty == '1']
NonPeine =  Crimes$homicide[Crimes$death_penalty == '0']
print(summary(OuiPeine))
print(summary(NonPeine))
ks.test(OuiPeine, NonPeine)

# Il faudrait superposer les deux histogrammes :
hist.default(OuiPeine, freq = FALSE)
hist.default(NonPeine, freq = FALSE)

# Test du chi2 d'independance entre la var homicide et blacks
chisq.test(Crimes$homicide, Crimes$white)
chisq.test(Crimes$homicide, Crimes$female)
# On veut savoir si le fait d'etre r??publicain est ind??pendant d'avoir la peine de mort
chisq.test(Crimes$female, Crimes$pauvrete)
# On affiche les scatterplots
vec = Crimes[,2:6]
pairs(vec)
# On superpose les histogrammes

set.seed(42)
p1 <- hist(Crimes$homicide)                     # centered at 4
p2 <- hist(Crimes$black)                     # centered at 6
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,100))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,100), add=T)  # second

#rob
############################################
## fonction pour l'acp centree reduite
############################################

myAcpRed <- function(data, graph = FALSE){
  
  
  #############################################
  ##Transformation des donnÃ©es et matrice de Variance
  #############################################
  
  dataMatrix <- as.matrix(data)
  cenRed <- function(x){y <- sqrt((1/nrow(dataMatrix))*sum((x - mean(x))^2))
  out <- (x - mean(x))/y
  return(out)}
  dataNew <- apply(dataMatrix, 2, cenRed)
  V <- (1/nrow(dataNew))*(t(dataNew)%*%dataNew)
  
  ###########################################
  ## recherche des composantes principales
  ###########################################
  
  diagoV <- eigen(V)
  
  composantePrincipal <- apply(diagoV[[2]],2, function(x){out <- dataNew%*%x
  return(out)})
  valeurPrincipal <- diagoV[[1]]
  inertieTotale <- sum(valeurPrincipal)
  Pourcentage_inertie <- valeurPrincipal/inertieTotale
  
  #########################################################################################################
  
  ###########################################
  ##Illustration avec les deux axes principaux (PC1 et PC2)
  ###########################################
  
  correlations <- matrix(c(sqrt(valeurPrincipal[1])*diagoV[[2]][,1], sqrt(valeurPrincipal[2])*diagoV[[2]][,2]), ncol = 2)
  
  if(graph){
    
    ############################################
    ##cercles des corrÃ©lations
    ############################################
    
    x11()
    u <- seq(0,2*pi,0.001)
    plot(cos(u), sin(u), type = "l", main = "cercle des corrÃ©lations pour le premier plan factoriel", xlab = "PC1", ylab = "PC2")
    polygon(c(0,0),c(-1,1),  lty = 3)
    polygon(c(-1,1), c(0,0), lty = 3)
    
    aux1 <- sapply(1:ncol(dataNew), function(int){arrows(0,0,correlations[int,1],correlations[int,2], angle = 15, length = 0.10)
      text(correlations[int, 1],correlations[int,2], labels = colnames(dataNew)[int])})
    
    ###############################################
    ##Projection sur le premier plan factoriel
    ###############################################
    
    x11()
    xx <- c(min(composantePrincipal[,1]), max(composantePrincipal[,1]))
    yy <- c(min(composantePrincipal[,2]), max(composantePrincipal[,2]))
    
    
    plot(composantePrincipal, type = "p", pch = 16, xlab  = "PC1", ylab = "PC2", 
         main = "Projections des individus sur le premier plan factoriel")
    polygon(c(0,0), c(yy[1],yy[2]), lty = 3)
    polygon(c(xx[1],xx[2]), c(0,0), lty = 3)
    
    aux2 <- sapply(1:nrow(dataNew), function(int){text(composantePrincipal[int,1], composantePrincipal[int,2]-0.15, labels = int)})
    
    
    x11()
    xx <- c(min(composantePrincipal[,1]), max(composantePrincipal[,1]))
    yy <- c(min(composantePrincipal[,2]), max(composantePrincipal[,2]))
    
    
    plot(composantePrincipal, type = "p", pch = 16, xlab  = "PC1", ylab = "PC2", main = "RÃ©sumÃ© sur le premier plan factoriel")
    polygon(c(0,0), c(yy[1],yy[2]), lty = 3)
    polygon(c(xx[1],xx[2]), c(0,0), lty = 3)
    
    u <- seq(0,2*pi,0.001)
    lines(cos(u), sin(u))
    
    
    aux2 <- sapply(1:nrow(dataNew), function(int){text(composantePrincipal[int,1], composantePrincipal[int,2]-0.15, labels = int)})
    aux1 <- sapply(1:ncol(dataNew), function(int){arrows(0,0,correlations[int,1],correlations[int,2], angle = 15, length = 0.10)
      text(correlations[int, 1],correlations[int,2], labels = colnames(dataNew)[int])})
    
    
    
    
    
    
  }
  
  return(list(composantePrincipal = composantePrincipal, valeurPrincipal = valeurPrincipal, 
              Pourcentage_inertie = Pourcentage_inertie, inertieTotale = inertieTotale, correlations = correlations))
}

#trouver les variables pivots
quant<-list(Crimes[1],Crimes[2],Crimes[4],Crimes[4])
acp_crimes<-prcomp(Crimes[,3:6,])
#acp_crimes<-prcomp(quant[,1:3,])
names(acp_crimes)
acp_crimes[2]
plot(acp_crimes)

#variable pivot
x1=-0.5171896255*Crimes$black +0.8558705549*Crimes$white+0.0005930161*Crimes$pauvrete+0.0003646115*Crimes$female
x2=0.5164854334*Crimes$white+0.8547422208*Crimes$black+0.0515578501*Crimes$female+0.0005666858*Crimes$pauvrete
reg3<-lm(Crimes$homicide~x1+x2+Crimes$revenu+Crimes$chomage+Crimes$death_penalty)
summary(reg3)

#acp2_crimes<-myAcpRed(Crimes[,3:8,],graph=T)
#plot(acp2_crimes)
