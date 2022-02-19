##SPON2 expression appears to correlate with tumor size in heatmaps - let's test this

##import dNF and ENCODE skin sample gene expression data 
source("../2016-11-17/dNF_ENCODE_GO_Mast.R")

##all gene sets pulled from Gene Ontology

##plot SPON2 expression against tumor size
SPON2 <- c("SPON2")
SPON2.df <- dplyr::filter(norm.df, Gene %in% SPON2)

tumorsize.df$Sample <- rownames(tumorsize.df)

SPON2.df <- dplyr::left_join(SPON2.df, tumorsize.df, by = "Sample")
SPON2.df <- dplyr::slice(SPON2.df, 1:33)

library(ggplot2)

SPON2.plot <- qplot(SPON2.df$Expression, SPON2.df$Size)
SPON2.plot + stat_smooth(method = "lm", formula = y ~ x, size = 1)

ggsave("SPON2.pdf")

##plot FCER1G against tumor size
FCER1G <- c("FCER1G")
FCER1G.df <- dplyr::filter(norm.df, Gene %in% FCER1G)

tumorsize.df$Sample <- rownames(tumorsize.df)

FCER1G.df <- dplyr::left_join(FCER1G.df, tumorsize.df, by = "Sample")
FCER1G.df <- dplyr::slice(FCER1G.df, 1:33)

library(ggplot2)

FCER1G.plot <- qplot(FCER1G.df$Expression, FCER1G.df$Size)
FCER1G.plot + stat_smooth(method = "lm", formula = y ~ x, size = 1)

ggsave("FCER1G.pdf")

##plot FCER1A against tumor size
FCER1A <- c("FCER1A")
FCER1A.df <- dplyr::filter(norm.df, Gene %in% FCER1A)

tumorsize.df$Sample <- rownames(tumorsize.df)

FCER1A.df <- dplyr::left_join(FCER1A.df, tumorsize.df, by = "Sample")
FCER1A.df <- dplyr::slice(FCER1A.df, 1:33)

library(ggplot2)

FCER1A.plot <- qplot(FCER1A.df$Expression, FCER1A.df$Size)
FCER1A.plot + stat_smooth(method = "lm", formula = y ~ x, size = 1)

ggsave("FCER1A.pdf")

##plot FGR against tumor size
FGR <- c("FGR")
FGR.df <- dplyr::filter(norm.df, Gene %in% FGR)

tumorsize.df$Sample <- rownames(tumorsize.df)

FGR.df <- dplyr::left_join(FGR.df, tumorsize.df, by = "Sample")
FGR.df <- dplyr::slice(FGR.df, 1:33)

library(ggplot2)

FGR.plot <- qplot(FGR.df$Expression, FGR.df$Size)
FGR.plot + stat_smooth(method = "lm", formula = y ~ x, size = 1)

ggsave("FGR.pdf")


##plot UNC13D against tumor size
UNC13D <- c("UNC13D")
UNC13D.df <- dplyr::filter(norm.df, Gene %in% UNC13D)

tumorsize.df$Sample <- rownames(tumorsize.df)

UNC13D.df <- dplyr::left_join(UNC13D.df, tumorsize.df, by = "Sample")
UNC13D.df <- dplyr::slice(UNC13D.df, 1:33)

library(ggplot2)

UNC13D.plot <- qplot(UNC13D.df$Expression, UNC13D.df$Size)
UNC13D.plot + stat_smooth(method = "lm", formula = y ~ x, size = 1)

ggsave("UNC130.pdf")

##plot KITLG against tumor size
KITLG <- c("KITLG")
KITLG.df <- dplyr::filter(norm.df, Gene %in% KITLG)

tumorsize.df$Sample <- rownames(tumorsize.df)

KITLG.df <- dplyr::left_join(KITLG.df, tumorsize.df, by = "Sample")
KITLG.df <- dplyr::slice(KITLG.df, 1:33)

library(ggplot2)

KITLG.plot <- qplot(KITLG.df$Expression, KITLG.df$Size)
KITLG.plot + stat_smooth(method = "lm", formula = y ~ x, size = 1)

ggsave("KITLG.pdf")

##plot C3 against tumor size
C3 <- c("C3")
C3.df <- dplyr::filter(norm.df, Gene %in% C3)

tumorsize.df$Sample <- rownames(tumorsize.df)

C3.df <- dplyr::left_join(C3.df, tumorsize.df, by = "Sample")
C3.df <- dplyr::slice(C3.df, 1:33)

library(ggplot2)

C3.plot <- qplot(C3.df$Expression, C3.df$Size)
C3.plot + stat_smooth(method = "lm", formula = y ~ x, size = 1)

ggsave("C3.pdf")


##plot FYN against tumor size
FYN <- c("FYN")
FYN.df <- dplyr::filter(norm.df, Gene %in% FYN)

tumorsize.df$Sample <- rownames(tumorsize.df)

FYN.df <- dplyr::left_join(FYN.df, tumorsize.df, by = "Sample")
FYN.df <- dplyr::slice(FYN.df, 1:33)

library(ggplot2)

FYN.plot <- qplot(FYN.df$Expression, FYN.df$Size)
FYN.plot + stat_smooth(method = "lm", formula = y ~ x, size = 1)

ggsave("FYN.pdf")


##plot NGF against tumor size
NGF <- c("NGF")
NGF.df <- dplyr::filter(norm.df, Gene %in% NGF)

tumorsize.df$Sample <- rownames(tumorsize.df)

NGF.df <- dplyr::left_join(NGF.df, tumorsize.df, by = "Sample")
NGF.df <- dplyr::slice(NGF.df, 1:33)

library(ggplot2)

NGF.plot <- qplot(NGF.df$Expression, NGF.df$Size)
NGF.plot + stat_smooth(method = "lm", formula = y ~ x, size = 1)

ggsave("NGF.pdf")


