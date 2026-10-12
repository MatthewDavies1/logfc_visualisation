
library(ggplot2)
library(dplyr)
library(ggrepel)

raw.buffer <- read.csv("inputR_buffer.txt", header = TRUE, sep = "\t")
raw.serum <- read.csv("inputR_serum.txt", header = TRUE, sep = "\t")

#######buffer 
input.buffer <- raw.buffer
input.buffer$colour <- ifelse(input.buffer$Adj..p.value < 0.05,
                              ifelse(input.buffer$log2FC >=1 , "#42a612",
                                     ifelse(input.buffer$log2FC <= -1 , "#a11414",
                                            ifelse(input.buffer$log2FC >=0.5 & input.buffer$log2FC <1, "#d38210",
                                                   ifelse(input.buffer$log2FC <= -0.5 & input.buffer$log2FC > -1, "#d38210",
                                                          "grey")))),
                              "grey")
input.buffer$size <- ifelse(input.buffer$Adj..p.value < 0.05,
                            ifelse(input.buffer$log2FC >=1 , 5,
                                   ifelse(input.buffer$log2FC <= -1, 5,
                                          ifelse(input.buffer$log2FC >= 0.5 & input.buffer$log2FC < 1, 5,
                                                 ifelse(input.buffer$log2FC <= -0.5 & input.buffer$log2FC > -1, 5,
                                                        2)))),
                            2)

lines_df <- data.frame(
  y = c(0, 1, -1, 0.5, -0.5),
  colour = c("black", "#848484", "#848484", "#848484", "#848484"),
  linetype = c("solid", "dashed", "dashed", "dashed", "dashed")
)

plot.buffer <- ggplot(input.buffer, aes(x = Orf, y = log2FC, colour = colour, size = size)) + 
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(2,5), guide = "none")+
  scale_y_continuous(limits = c(-2.5,2.5), breaks = seq(-2,2,1)) + 
  geom_hline(aes(yintercept = y, color = colour, linetype = linetype), linewidth = c(1, 0.1, 0.1, 0.1, 0.1), data = lines_df) +
  scale_color_identity() + 
  scale_linetype_identity() +
  ylab("log2FC") +
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 16)) +
  geom_text_repel(data = input.buffer, box.padding = 0.7, 
                  max.overlaps = Inf, size=6,
                  aes(label=annotate, colour = "black", fontface = 'italic'))+
  geom_text_repel(data = input.buffer, box.padding = 0.7, 
                  max.overlaps = Inf, size=6,
                  aes(label=annotate2, colour = "black"))

plot.buffer
ggsave("plot_logfc_buffer.png", plot = plot.buffer, width = 16, height = 11, units = "in", dpi = 500)
ggsave("plot_logfc_buffer.svg", plot = plot.buffer, width = 16, height = 11, units = "in", dpi = 500)

#######serum
input.serum <- raw.serum
input.serum$colour <- ifelse(input.serum$Adj.P.Value < 0.05,
                              ifelse(input.serum$Log.2.FC >=1 , "#42a612",
                                     ifelse(input.serum$Log.2.FC <= -1 , "#a11414",
                                            ifelse(input.serum$Log.2.FC >=0.5 & input.serum$Log.2.FC <1, "#d38210",
                                                   ifelse(input.serum$Log.2.FC <= -0.5 & input.serum$Log.2.FC > -1, "#d38210",
                                                          "grey")))),
                              "grey")
input.serum$size <- ifelse(input.serum$Adj.P.Value < 0.05,
                            ifelse(input.serum$Log.2.FC >=1 , 5,
                                   ifelse(input.serum$Log.2.FC <= -1, 5,
                                          ifelse(input.serum$Log.2.FC >= 0.5 & input.serum$Log.2.FC < 1, 5,
                                                 ifelse(input.serum$Log.2.FC <= -0.5 & input.serum$Log.2.FC > -1, 5,
                                                        2)))),
                            2)


plot.serum <- ggplot(input.serum, aes(x = X.ORF, y = Log.2.FC, colour = colour, size = size)) + 
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(2,5), guide = "none")+
  scale_y_continuous(limits = c(-2.5,2.5), breaks = seq(-2,2,1)) + 
  geom_hline(aes(yintercept = y, color = colour, linetype = linetype), linewidth = c(1, 0.5, 0.5, 0.5, 0.5), data = lines_df) +
  scale_color_identity() + 
  scale_linetype_identity() +
  ylab("log2FC") +
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 16))+
  geom_text_repel(data = input.serum, box.padding = 0.5, 
                  max.overlaps = Inf, size=6,
                  aes(label=annotate, colour = "black", fontface = 'italic'))+
  geom_text_repel(data = input.serum, box.padding = 0.7, 
                  max.overlaps = Inf, size=6,
                  aes(label=annotate2, colour = "black"))

plot.serum
ggsave("plot_logfc_serum.png", plot = plot.serum, width = 16, height = 11, units = "in", dpi = 500)
ggsave("plot_logfc_serum.svg", plot = plot.serum, width = 16, height = 11, units = "in", dpi = 500)

