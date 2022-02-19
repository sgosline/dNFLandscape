library(synapseClient)
library(deconstructSigs)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
synapseLogin()

mut <- fread(synGet("syn5713423")@filePath)
mut2 <- fread(synGet("syn6087182")@filePath)

sigs <- mut.to.sigs.input(mut.ref = mut,
                          sample.id = "Sample_ID",
                          chr = "Chromosome",
                          pos = "Start_Position",
                          ref = "Reference_Allele",
                          alt = "Variant_Allele")

sigs2 <- mut.to.sigs.input(mut.ref = mut2,
                          sample.id = "Sample",
                          chr = "Chr",
                          pos = "Start",
                          ref = "Ref",
                          alt = "Alt")

#for(i in rownames(sigs)){
#print(i)
#x<-whichSignatures(tumor.ref = sigs,
#                signatures.ref = signatures.cosmic,
#                sample.id = i,
#                contexts.needed = TRUE,
#                tri.counts.method = "default")
#png(paste(i,"cosmicsig_bar.png", sep = ""))
#plotSignatures(x, sub = i)
#dev.off()
#
#
#png(paste(i,"cosmicsig_pie.png", sep = ""))
#makePie(x)
#dev.off()
#}

#for(i in rownames(sigs2)){
#  print(i)
#  x<-whichSignatures(tumor.ref = sigs2,
#                     signatures.ref = signatures.cosmic,
#                     sample.id = i,
#                     contexts.needed = TRUE,
#                     tri.counts.method = "default")
#  png(paste(i,"syn6047149_cosmicsig_bar.png", sep = ""))
#  plotSignatures(x, sub = i)
#  dev.off()
  
#  png(paste(i,"syn6047149_cosmicsig_pie.png", sep = ""))
#  makePie(x)
#  dev.off()
#}

#x<-whichSignatures(tumor.ref = sigs,
#                   signatures.ref = signatures.nature2013,
#                   sample.id = "patient_4_tissue_0005",
#                   contexts.needed = TRUE,
#                   tri.counts.method = "default")
#plotSignatures(x)
#makePie(x)

foo<-lapply(rownames(sigs), function(i){
  x<-whichSignatures(tumor.ref = sigs,
                     signatures.ref = signatures.cosmic,
                     sample.id = i,
                     contexts.needed = TRUE,
                     tri.counts.method = "default")
  x$weights
})
status <- distinct(select(mut, Mutation_Status, Sample_ID))
names(status) <- c("mut_stat", ".id")
names(foo) <- rownames(sigs)
foo <- ldply(foo)
foo <- gather(foo, "Signature", "weight", 2:31)
foo$patient <- gsub("_tissue_....","",foo$.id)
foo <- left_join(foo, status)

foo.som <- foo %>% filter(mut_stat == "Somatic")
bar <- foo.som %>% group_by(Signature) %>% summarise(sum(weight))
foo.som <- left_join(foo.som, bar)
foo.som$sum <- foo.som$`sum(weight)`


foo.germ <- foo %>% filter(mut_stat == "Germline")
bar <- foo.germ %>% group_by(Signature) %>% summarise(sum(weight))
foo.germ <- left_join(foo.germ, bar)
foo.germ$sum <- foo.germ$`sum(weight)`


ggplot(data = foo.germ %>% filter(weight > 0.0), 
       aes(x = Signature %>% reorder(-sum), y = weight, fill = .id)) + 
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x=element_text(angle=60,hjust=1)) +
  labs(x= "Signature", y = "additive weight", title = "Germline Variant Signatures")
ggsave("additive_germline_sig_syn5713423.png")

ggplot(data = foo.germ %>% filter(weight > 0.0), 
       aes(x = .id, y = weight, fill = Signature)) + 
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x=element_text(angle=60,hjust=1)) +
  labs(x= "Signature", y = "weight", title = "Germline Variant Signatures")
ggsave("byPatient_germline_sig_syn5713423.png")

ggplot(data = foo.som %>% filter(weight > 0.0), 
       aes(x = Signature %>% reorder(-sum), y = weight, fill = .id)) + 
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x=element_text(angle=60,hjust=1)) +
  labs(x= "Signature", y = "additive weight", title = "Somatic Variant Signatures")
ggsave("additive_somatic_sig_syn5713423.png")

ggplot(data = foo.som %>% filter(weight > 0.0), 
       aes(x = .id, y = weight, fill = Signature)) + 
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x=element_text(angle=60,hjust=1)) +
  labs(x= "Signature", y = "weight", title = "Somatic Variant Signatures")
ggsave("byPatient_somatic_sig_syn5713423.png")



foo2<-lapply(rownames(sigs2), function(i){
  x<-whichSignatures(tumor.ref = sigs2,
                     signatures.ref = signatures.cosmic,
                     sample.id = i,
                     contexts.needed = TRUE,
                     tri.counts.method = "default")
  x$weights
})

names(foo2) <- rownames(sigs2)
foo2 <- ldply(foo2)
foo2 <- gather(foo2, "Signature", "weight", 2:31)

ggplot(data = foo2 %>% filter(weight > 0), aes(x = .id, y = weight, fill = Signature)) + 
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x=element_text(angle=60,hjust=1)) +
  labs(x= "Signature", y = "additive weight", title = "Somatic Variant Signatures")
ggsave("byPatient_somatic_sig_syn6087182.png")

ggplot(data = foo2 %>% filter(weight > 0), aes(x = Signature, y = weight, fill = .id)) + 
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x=element_text(angle=60,hjust=1)) +
  labs(x= "Signature", y = "additive weight", title = "Somatic Variant Signatures")
ggsave("additive_somatic_sig_syn6087182.png")
