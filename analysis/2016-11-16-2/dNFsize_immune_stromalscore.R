##library(synapseClient)
##synapseLogin()
##purity.scores <- synGet("syn6040997", load = TRUE)

purity.scores <- read.delim("estimateScoresWithOtherTumorData.txt")

purity.scores$Patient <- as.factor(purity.scores$Patient) 

library(ggplot2)

##plot tumor size vs Stromal score without color coding, linear fit and loess fit
size.Stromal <- qplot(purity.scores$TumorSize, purity.scores$StromalScore)
size.Stromal + stat_smooth(method = "lm", formula = y ~ x, size = 1)
size.Stromal + stat_smooth(method = "loess", formula = y ~ x, size = 1)

##plot tumor size vs Stromal score with color coding by patient number, linear fit only
size.Stromal <- qplot(purity.scores$TumorSize, purity.scores$StromalScore)
size.Stromal + stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  geom_point(aes(color=purity.scores$Patient), size = 3) +
  scale_color_brewer(palette = "Paired")

##plot tumor size vs Immune score without color coding, linear fit and loess fit
size.Immune <- qplot(purity.scores$TumorSize, purity.scores$ImmuneScore)
size.Immune + stat_smooth(method = "lm", formula = y ~ x, size = 1)
size.Immune + stat_smooth(method = "loess", formula = y ~ x, size = 1)

##plot tumor size vs Immune score with color coding by patient number, linear fit only
size.Immune <- qplot(purity.scores$TumorSize, purity.scores$ImmuneScore)
size.Immune + stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  geom_point(aes(color=purity.scores$Patient), size = 3) +
  scale_color_brewer(palette = "Paired")

##plot tumor size vs Purity score without color coding, linear fit and loess fit
size.Purity <- qplot(purity.scores$TumorSize, purity.scores$TumorPurity)
size.Purity + stat_smooth(method = "lm", formula = y ~ x, size = 1)
size.Purity + stat_smooth(method = "loess", formula = y ~ x, size = 1)

##plot tumor size vs Purity score with color coding by patient number, linear fit only
size.Purity <- qplot(purity.scores$TumorSize, purity.scores$TumorPurity)
size.Purity + stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  geom_point(aes(color=purity.scores$Patient), size = 3) +
  scale_color_brewer(palette = "Paired")

##plot tumor size vs ESTIMATE score without color coding, linear fit and loess fit
size.ESTIMATEscore <- qplot(purity.scores$TumorSize, purity.scores$ESTIMATEScore)
size.ESTIMATEscore + stat_smooth(method = "lm", formula = y ~ x, size = 1)
size.ESTIMATEscore + stat_smooth(method = "loess", formula = y ~ x, size = 1)

##plot tumor size vs ESTIMATE score with color coding by patient number, linear fit only
size.ESTIMATEscore <- qplot(purity.scores$TumorSize, purity.scores$ESTIMATEScore)
size.ESTIMATEscore + stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  geom_point(aes(color=purity.scores$Patient), size = 3) +
  scale_color_brewer(palette = "Paired")



