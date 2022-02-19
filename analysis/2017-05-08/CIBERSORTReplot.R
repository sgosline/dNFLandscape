library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(synapseClient)
library(RColorBrewer)
synapseLogin()

this.file <- "https://github.com/allaway/dNFLandscape/blob/master/analysis/2017-05-08/CIBERSORTReplot.R"

cib <- read.csv(synGet("syn5809355")@filePath, header = TRUE, stringsAsFactors = FALSE)
cib[1,1] <- "3096-PBK-0001"

rna.id<-read.table(synGet("syn6156140")@filePath, sep = ",", header = TRUE)
pt.id<-synTableQuery("SELECT 'RNASeq (Cufflinks)', sampleIdentifier FROM syn5556216")@values %>% filter(!is.na(`RNASeq (Cufflinks)`))
colnames(pt.id) <- c("SynapseID", "patientID")
ids <- full_join(rna.id,pt.id)

cib2 <- gather(cib, "Cell.Type", "CIBERSORT.Data", 2:23)

cib2$Cell.Type <- gsub("\\.", " ", cib2$Cell.Type)
cib2$Cell.Type <- gsub("T cells regulatory  Tregs ", "T cells regulatory Tregs", cib2$Cell.Type)
cib2$Cell.Type <- as.factor(cib2$Cell.Type)

map <- read.table(synGet("syn9766992")@filePath, header = TRUE, sep = "\t")
colnames(map) <- c("Cell.Type", "Cell.Plot.Name")
#cib2$Cell.Type <- factor(cib2$Cell.Type, levels = levels(reorder(cib2$Cell.Type , -cib2$CIBERSORT.Data, FUN = median)))
#cib2$median <- ave(cib2$CIBERSORT.Data, as.factor(cib2$Cell.Type), FUN=median)

cib2 <- full_join(cib2, map)
##sorting by mean 
cib2$Cell.Plot.Name <- factor(cib2$Cell.Plot.Name, levels = levels(reorder(cib2$Cell.Plot.Name , cib2$CIBERSORT.Data, FUN = mean)))
cib2$mean<- ave(cib2$CIBERSORT.Data, as.factor(cib2$Cell.Plot.Name), FUN=mean)

ggplot(cib2, aes(x=Cell.Plot.Name, y=CIBERSORT.Data, fill = mean)) +
  geom_boxplot(color = "black") + 
  scale_fill_gradient(low = "#004E98", high = "#FE9920") +
  theme(legend.position="none",
        text = element_text(size=15),
        axis.text.x=element_text(angle=60,hjust=1),
        plot.margin=unit(c(1,1,1,1), "cm")) +
  labs(x = "Immune Cell Type", y = "CIBERSORT Output (Fraction)") +
  coord_flip()

ggsave("cibersort_review_figure.png")
synStore(File("cibersort_review_figure.png", parentId = "syn5809348"), used = c("syn5556216", "syn6156140", "syn5809355", "syn9766992"), executed = this.file)


