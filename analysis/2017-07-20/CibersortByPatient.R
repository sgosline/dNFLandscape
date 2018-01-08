source("../../bin/wgsAnalysis.R")
source("../../analysis/2016-09-20/geneSampleMatrix.R")
source("../../bin/ssGSEAandGSVA_perPatient.R")
library(ggplot2)
library(ggbeeswarm)

cib.scores<-read.table(synGet('syn5809355')@filePath,sep=',',header=T)
cib.scores$Input.Sample<-as.character(cib.scores$Input.Sample)
cib.scores[1,1] <- "3096-PBK-0001"
colnames(cib.scores)[1] <- "Sample.Id"

ids<-read.table(synGet("syn6156140")@filePath, sep = ",", header = T)

colnames(ids)[4] <- "RNASeq (Cufflinks)"
moreids<-synTableQuery("SELECT * FROM syn5556216")@values %>% filter(!is.na(RNASeq)) %>% 
  left_join(ids) %>% 
  dplyr::select(sampleIdentifier, Sample.Id) %>% 
  arrange(Sample.Id)

cib.scores<-cib.scores %>% left_join(moreids)
rownames(cib.scores) <- cib.scores$sampleIdentifier

cib.scores$patient <- as.character(t(as.data.frame(strsplit(cib.scores$sampleIdentifier, split = "tumor"))[1,]))

cib.scores <- cib.scores %>% group_by(patient) %>% 
  mutate("var.mast" = var(Mast.cells.resting), "var.mac" = var(Macrophages.M2))

ggplot(data = cib.scores) +
  geom_boxplot(aes(x= patient, y= Mast.cells.resting, fill = patient), color = "black")+
  geom_beeswarm(aes(x= patient, y= Mast.cells.resting))+
  geom_text(aes(x = patient, y = 0.7, label = signif(var.mast, 2)), angle = 45)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(limit = c(0, 0.8))
ggsave("MastCellsResting_Bypatient.png", height = 5, width = 5)


ggplot(data = cib.scores) +
  geom_boxplot(aes(x= patient, y= Macrophages.M2, fill = patient), color = "black")+
  geom_beeswarm(aes(x= patient, y= Macrophages.M2))+
  geom_text(aes(x = patient, y = 0.7, label = signif(var.mac, 2)), angle = 45)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_y_continuous(limit = c(0, 0.8))
ggsave("M2_Mac_Bypatient.png", height = 5, width = 5)


hallmarks <- as.data.frame(t(hallmark.ssGSEA))
hallmarks$tumor <- rownames(hallmarks)
hallmarks$patient <- t(as.data.frame(strsplit(hallmarks$tumor, split = "tumor")))[,1]

hallmarks <- hallmarks %>% group_by(patient) %>% 
  mutate("var.ifna" = var(HALLMARK_INTERFERON_ALPHA_RESPONSE), 
         "var.ifng" = var(HALLMARK_INTERFERON_GAMMA_RESPONSE),
         "var.comp" = var(HALLMARK_COMPLEMENT))

ggplot(data = hallmarks) +
  geom_boxplot(aes(x= patient, y= HALLMARK_INTERFERON_ALPHA_RESPONSE, fill = patient), color = "black")+
  geom_beeswarm(aes(x= patient, y= HALLMARK_INTERFERON_ALPHA_RESPONSE))+
  geom_text(aes(x = patient, y = -0.5, label = signif(var.ifna, 2)), angle = 45)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(limit = c(-0.6, 0.6))
ggsave("HALLMARK_IFNAlpha_perpatient.png", height = 4, width = 5)

ggplot(data = hallmarks) +
  geom_boxplot(aes(x= patient, y= HALLMARK_INTERFERON_GAMMA_RESPONSE, fill = patient), color = "black")+
  geom_beeswarm(aes(x= patient, y= HALLMARK_INTERFERON_GAMMA_RESPONSE))+
  geom_text(aes(x = patient, y = -0.5, label = signif(var.ifng, 2)), angle = 45)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(limit = c(-0.6, 0.6))
ggsave("HALLMARK_IFNGamma_perpatient.png", height = 4, width = 5)

ggplot(data = hallmarks) +
  geom_boxplot(aes(x= patient, y= HALLMARK_COMPLEMENT, fill = patient), color = "black")+
  geom_beeswarm(aes(x= patient, y= HALLMARK_COMPLEMENT))+
  geom_text(aes(x = patient, y = -0.5, label = signif(var.comp, 2)), angle = 45)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limit = c(-0.6, 0.6))
ggsave("HALLMARK_Complement_perpatient.png", height = 4, width = 5)


all.data <- full_join(hallmarks,cib.scores)

cor(all.data$HALLMARK_INTERFERON_ALPHA_RESPONSE, all.data$Macrophages.M2, method = "spearman")
cor(all.data$HALLMARK_INTERFERON_GAMMA_RESPONSE, all.data$Mast.cells.resting, method = "spearman")
cor(all.data$Mast.cells.resting, all.data$Macrophages.M2, method = "spearman")

ggplot(all.data) +
  geom_point(aes(x=HALLMARK_COMPLEMENT, y=Mast.cells.resting))


ihc<-read.table("../../data/Carroll_cNF_IHC.txt", sep = "\t", header = T) 
ihc$CD117 <- gsub("\\%", "", ihc$CD117)
ihc$Iba1 <- gsub("\\%", "", ihc$Iba1)
ihc$CD117 <- as.numeric(ihc$CD117)
ihc$Iba1 <- as.numeric(ihc$Iba1)
ihc$patient <- gsub("-\\d+-.+", "", ihc$X)
ihc$patient <- gsub("^0", "", ihc$patient)
ihc$patient <- gsub("^", "patient", ihc$patient)

mast.cib <- ungroup(all.data) %>% dplyr::select(patient, Mast.cells.resting, Macrophages.M2)

mast <- ihc %>% dplyr::select(patient, CD117) %>%
  inner_join(mast.gsva) %>% 
  distinct() %>% 
  filter(!is.na(CD117)) %>%
  mutate()


ggplot(data = ihc) +
  geom_boxplot(aes(x= patient, y= CD117, fill = patient), color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("MastCells_IHC.png", width = 5, height = 4)

ggplot(data = mast.cib %>% filter(patient %in% unique(ihc$patient))) +
  geom_boxplot(aes(x= patient, y= Mast.cells.resting, fill = patient), color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("MastCells_ssGSEA.png", width = 5, height = 4)

ggplot(data = ihc) +
  geom_boxplot(aes(x= patient, y= Iba1, fill = patient), color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Macrophages_IHC.png", width = 5, height = 4)

ggplot(data = mast.cib %>% filter(patient %in% unique(ihc$patient))) +
  geom_boxplot(aes(x= patient, y= Macrophages.M2, fill = patient), color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Macrophages_ssGSEA.png", width = 5, height = 4)
