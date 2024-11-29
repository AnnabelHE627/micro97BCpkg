###This is the code for R function by Sheng Xu

usethis::use_package(package="beeswarm",type="Imports")
usethis::use_package(package="vegan",type="Imports")
usethis::use_package(package="grid",type="Imports")
usethis::use_package(package="gridExtra",type="Imports")
usethis::use_package(package="ggplot2",type="Imports")
usethis::use_package(package="ggpubr",type="Imports")
usethis::use_package(package="dplyr",type="Imports")
usethis::use_package(package="reshape2",type="Imports")


#' Beeplot
#'
#' @param overweight data for overweight group
#' @param control data for control group
#' @param species microbiom species
#' @param group group variables
#' @param data dataframe
#' @param main2 main title
#'
#' @param xlab xlab
#' @param ylab ylab
#'
#' @return plot
#' @export
#' @importFrom beeswarm beeswarm
#'
#'
beeplot = function(overweight,control,species, group, data,main="",xlab="",ylab="")
{
  wilcoxTest <- wilcox.test(overweight,control,paired = T,alternative = "two.sided")$p.value

  if(wilcoxTest<0.001){
    wilcoxTest=signif(wilcoxTest,4)
    wilcoxTest=format(wilcoxTest,scientific = TRUE)
  }else{
    wilcoxTest=round(wilcoxTest,3)
  }
  ySeg=max(data[,species])*0.9
  boxplot(data[,species]~data[,group],data, main= main,
          cex.main=1.5,cex.lab=1.3,cex.axis=1.2,outline=T,border="black",ylab = ylab,xlab = "",names=c("Controls", "Cases"))#xlab=xlab,ylab=ylab)#,)x=group
  beeswarm(data[,species]~data[,group],data, col = c("blue","red"),lwd=0.2,
           pch = 16, add = TRUE, corral="wrap")
  segments(1,ySeg, 2,ySeg);segments(1,ySeg, 1,ySeg*0.96);segments(2,ySeg, 2,ySeg*0.96)
  text(1.5,ySeg*1.05,labels=paste("p=",wilcoxTest,sep=""))
}




beeplot_unpaired = function(species, group, data)
{
  wilcoxTest <- wilcox.test(data[,species]~data[,group],data)$p.value

  if(wilcoxTest<0.001){
    wilcoxTest=signif(wilcoxTest,4)
    wilcoxTest=format(wilcoxTest,scientific = TRUE)
  }else{
    wilcoxTest=round(wilcoxTest,3)
  }
  ySeg=max(data[,species])*0.9
  boxplot(data[,species]~data[,group],data, main= "",
          cex.main=1.5,cex.lab=1.3,cex.axis=1.2,outline=T,border="black",xlab = group,ylab = species)
  beeswarm(data[,species]~data[,group],data, col = c("blue","red"),lwd=0.2,
           pch = 16, add = TRUE, corral="wrap")
  segments(1,ySeg, 2,ySeg);segments(1,ySeg, 1,ySeg*0.96);segments(2,ySeg, 2,ySeg*0.96)
  text(1.5,ySeg*1.05,labels=paste("p=",wilcoxTest,sep=""))
}

#' Shannon_diversity
#'
#' @param metaplan_species_all data for microbiome species
#'
#' @return plot
#' @export
#'
#' @importFrom vegan diversity
#'
Shannon_diversity<-function(metaplan_species_all){
  Shannon_diversity <- metaphlan_species_all

  Shannon_diversity <- as.data.frame(vegan::diversity(Shannon_diversity, index = "shannon", MARGIN = 2, base = exp(1)))

  Shannon_diversity$ID <- row.names(Shannon_diversity)

  Shannon_diversity$Number <- Shannon_diversity$`vegan::diversity(Shannon_diversity, index = "shannon", MARGIN = 2, base = exp(1))`


  Shannon_diversity$Overweight <-  metaphlan_species[match(rownames(Shannon_diversity),rownames(metaphlan_species)),"Overweight"]

  Shannon_diversity$sex <-  metaphlan_species[match(rownames(Shannon_diversity),rownames(metaphlan_species)),"sex"]
  Shannon_diversity$subclass <-  metaphlan_species[match(rownames(Shannon_diversity),rownames(metaphlan_species)),"subclass"]
  Shannon_diversity <- Shannon_diversity[order(Shannon_diversity$subclass),]
  Shannon_diversity <- Shannon_diversity[-which(is.na(Shannon_diversity$Overweight)),]

  overweight <- subset(Shannon_diversity,  Overweight == "1", Number,drop = TRUE)
  control <- subset(Shannon_diversity,  Overweight == "0", Number,drop = TRUE)

  wilcox.test(overweight,control,paired = T,alternative = "two.sided")

  beeplot(overweight,control,"Number","Overweight",Shannon_diversity)

}

#' Title
#'
#' @param metaplan_species_all
#'
#' @return plot
#' @export
#'
#' @importFrom vegan diversity
simpson_diversity<-function(metaplan_species_all){
  simpson_diversity <- metaphlan_species_all

  simpson_diversity <- as.data.frame(vegan::diversity(simpson_diversity, index = "invsimpson", MARGIN = 2, base = exp(1)))

  simpson_diversity$ID <- row.names(simpson_diversity)

  simpson_diversity$Index <- simpson_diversity$`vegan::diversity(simpson_diversity, index = "invsimpson", MARGIN = 2, base = exp(1))`


  simpson_diversity$Overweight <-  metaphlan_species[match(rownames(simpson_diversity),rownames(metaphlan_species)),"Overweight"]

  simpson_diversity$sex <-  metaphlan_species[match(rownames(simpson_diversity),rownames(metaphlan_species)),"sex"]
  simpson_diversity$subclass <-  metaphlan_species[match(rownames(simpson_diversity),rownames(metaphlan_species)),"subclass"]
  simpson_diversity <- simpson_diversity[order(simpson_diversity$subclass),]
  simpson_diversity <- simpson_diversity[-which(is.na(simpson_diversity$Overweight)),]

  overweight_ivSimpson <- subset(simpson_diversity,  Overweight == "1", Index,drop = TRUE)
  control_ivSimpson <- subset(simpson_diversity,  Overweight == "0", Index,drop = TRUE)

  wilcox.test(overweight_ivSimpson,control_ivSimpson,paired = T,alternative = "two.sided")

  simpson_diversity$Group<-ifelse(simpson_diversity$Overweight==1,"Overweight","Non-overweight")
  beeplot(overweight_ivSimpson,control_ivSimpson,"Index","Group",simpson_diversity, main="InvSimpson")
}

#' Title
#'
#' @param metaplan_species_all
#'
#' @return plot
#' @export
#'
#' @importFrom reshape2 melt
Pielou_evenness<-function(metaplan_species_all){

  Pielou_evenness <- metaphlan_species_all

  H <- diversity(t(metaphlan_species_all))
  S <- specnumber(t(metaphlan_species_all))
  Pielou_evenness <- as.data.frame(H/log(S))


  Pielou_evenness$Overweight <-  metaphlan_species[match(rownames(Pielou_evenness),rownames(metaphlan_species)),"Overweight"]
  Pielou_evenness$subclass <-  metaphlan_species[match(rownames(Pielou_evenness),rownames(metaphlan_species)),"subclass"]
  Pielou_evenness$index <- Pielou_evenness$`H/log(S)`
  Pielou_evenness <- Pielou_evenness[order(Pielou_evenness$subclass),]
  Pielou_evenness <- Pielou_evenness[-which(is.na(Pielou_evenness$Overweight)),]

  overweight_Pielou <- subset(Pielou_evenness,  Overweight == "1", index,drop = TRUE)
  control_Pielou <- subset(Pielou_evenness,  Overweight == "0", index,drop = TRUE)

  wilcox.test(overweight_Pielou,control_Pielou ,paired = T,alternative = "two.sided")


  Pielou_evenness$Group<-ifelse(Pielou_evenness$Overweight==1,"Overweight","Non-overweight")
  beeplot(overweight_Pielou,control_Pielou ,"index","Overweight",Pielou_evenness, main="Pielou evenness")
}
