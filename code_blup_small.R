datafile<-read.table('https://raw.githubusercontent.com/fragomen/lesson6.1/main/lesson6/datafile_effects_s.tsv',header=T)
snps<-read.table('https://raw.githubusercontent.com/fragomen/lesson6.1/main/lesson6/snp_file_effects_s.tsv',header=T)


dim(snps)



library(pedigree)
h2<-0.3
pedigree1 <- datafile %>% select(id,father,mother,pheno,sex,farm) %>%
  dplyr::mutate(farm=as.factor(farm),sex=as.factor(sex)) %>%  
  dplyr::rename(ID=id) %>% 
  as.data.frame() 
nanim<-dim(pedigree1)[1]
dim(datafile)
dim(snps)
#makeA(pedigree1,which = rep(TRUE))


blup_1<-blup(pheno~1,pedigree1,alpha=1/h2-1) #alpha sigma_e/sigma_a
# WE HAVE TO REMOVE FIXED EFFECT SOLUTIONS
ebv.temp<-(blup_1[-1,]) #MEAN

datafile2 <- datafile %>%  dplyr::mutate(ebv=ebv.temp)


ebv.temp1<-ebv.temp %>%  as_tibble() %>%
  dplyr::mutate(id=names_blup) %>%
  dplyr::rename(ebv=value)

datafile2 %>% summarize(acc_blup=cor(gv,ebv))

########now the GBLUP
rownames(snps) <- pedigree1$ID[2001:4000]
snp2 <- as.matrix(snps)
rownames(snp2)
head(pedigree1)
gblup_1<-gblup(pheno~1,data=pedigree1,M=snp2,lambda=1/h2-1) #alpha sigma_e/sigma_a



gebv.temp<-(gblup_1[-1,]) %>%  as_tibble() %>%
  dplyr::rename(gebv=value)

datafile3 <- datafile2 %>%  as_tibble() %>% dplyr::filter(gen>5) %>% 
  dplyr::mutate(gebv=gebv.temp)


datafile3 %>% summarize(acc_blup=cor(gv,ebv),acc_gblup=cor(gv,gebv),
                        cor_ebv_gebv=cor(ebv,gebv))




######################
##### Reduced Data  ##
######################

datafile_reduced <- datafile %>% dplyr::mutate(pheno2=ifelse(gen>9,NA,pheno))
pedigree2 <- datafile_reduced %>% select(id,father,mother,pheno2)%>%
  dplyr::rename(ID=id) %>% 
  as.data.frame()


blup_2<-blup(pheno2~1,pedigree2,alpha=1/h2-1) #alpha sigma_e/sigma_a
# WE HAVE TO REMOVE FIXED EFFECT SOLUTIONS
ebv.temp<-(blup_2[-1,]) #MEAN

datafile2_red <- datafile_reduced %>%  dplyr::mutate(ebv=ebv.temp)


ebv.temp1<-ebv.temp %>%  as_tibble() %>%
  dplyr::mutate(id=names_blup) %>%
  dplyr::rename(ebv=value)

datafile2_red %>% summarize(acc_blup=cor(gv,ebv))

########now the GBLUP
rownames(snps) <- pedigree2$ID[2001:4000]
snp2 <- as.matrix(snps)
rownames(snp2)

gblup_2<-gblup(pheno2~1,data=pedigree2,M=snp2,lambda=1/h2-1) #alpha sigma_e/sigma_a



gebv.temp<-(gblup_2[-1,]) %>%  as_tibble() %>%
  dplyr::rename(gebv=value)

datafile3_reduced <- datafile2_red %>%  as_tibble() %>%
  dplyr::filter(gen>5) %>% 
  dplyr::mutate(gebv=gebv.temp)


datafile3_reduced %>% summarize(acc_blup=cor(gv,ebv),acc_gblup=cor(gv,gebv),
                                cor_ebv_gebv=cor(ebv,gebv))


#calculating accuracy for young animals only

datafile3_reduced %>%   dplyr::filter(gen==10) %>% 
  summarize(acc_blup=cor(gv,ebv),acc_gblup=cor(gv,gebv),
            cor_ebv_gebv=cor(ebv,gebv))

  #homework: test fix effects in the model :)