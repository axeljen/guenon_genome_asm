rm(list = ls())
library(tidyverse)

# get syntenymap, genomes and calls from command line
# define args as input argument vector
args <- commandArgs(trailingOnly = TRUE)

# chromcords as first as syntney map as second argument
synmap_file <- args[1]
genomes_file <- args[2]
calls_file <- args[3]
outpref <- args[4]

# syntenymap
synmap <- read_tsv(synmap_file)

# genomes 
genomes <- read_tsv(genomes_file)

# fissions and fusions
calls <- read_tsv(calls_file)

fissions <- calls %>%
  filter(type == "fission") %>%
  mutate(target_min_cum = target_pos_1_min + genomes$cumulative_start[match(target_chrom_1,genomes$chrom)],
         target_max_cum = target_pos_1_max + genomes$cumulative_start[match(target_chrom_1,genomes$chrom)]) 

fusions <- calls %>%
  filter(type == "fusion") %>%
  mutate(target_pos_1_cum = target_pos_1_min + genomes$cumulative_start[match(target_chrom_1,genomes$chrom)],
         target_pos_2_cum = target_pos_2_min + genomes$cumulative_start[match(target_chrom_2,genomes$chrom)])


# add ypositions to genomes
genomes <- genomes %>%
  mutate(ypos = ifelse(genome == "target", 10,0))

# and synmap
synmap <- synmap %>%
  mutate(target_y = genomes$ypos[match(tchrom,genomes$chrom)],
         query_y = genomes$ypos[match(qchrom,genomes$chrom)]) 

# and fissions/fusions
fissions <- fissions %>%
  mutate(ypos = genomes$ypos[match(target_chrom_1,genomes$chrom)])

fusions <- fusions %>%
  mutate(ypos = genomes$ypos[match(target_chrom_1,genomes$chrom)])

# make plygonized synmap positions
synmap_long <- synmap %>%
  mutate(x1 = tstart_cumulative, x2 = tend_cumulative,
         x3 = tend_cumulative, x4 = qend_cumulative,
         x5 = qend_cumulative,
         x6 = qstart_cumulative, x7 = qstart_cumulative,
         x8 = tstart_cumulative, x9 = tstart_cumulative,
         y1 = target_y - 1, y2 = target_y -1,
         y3 = target_y - 3, y4 = query_y + 3, y5 = query_y + 1,
         y6 = query_y + 1, y7 = query_y + 3, y8 = target_y - 3, y9 = target_y - 1) %>%
  mutate(alignment = as.character(rownames(.))) %>%
  pivot_longer(cols = c(x1, x2, x3, x4, x5, x6, x7, x8, x9), names_to = "x_col", values_to = "x") %>%
  pivot_longer(cols = c(y1, y2, y3, y4, y5, y6, y7, y8, y9), names_to = "y_col", values_to = "y") %>%
  # Ensure that the 'x' and 'y' columns align properly
  group_by(alignment) %>%
  filter(str_extract(x_col, "\\d") == str_extract(y_col, "\\d")) %>%
  ungroup()


# plot the genomes
p <- ggplot() +
  geom_rect(data = genomes, aes(xmin = cumulative_start, xmax = cumulative_start + length,
                ymin = ypos -1, ymax = ypos + 1), color = "black",
            fill = "white") +
  # add polygons for the alignments
  #geom_segment(data = synmap, aes(x = tstart_cumulative, xend = qstart_cumulative,
                               #y = target_y - 1, yend = query_y + 1))
  geom_polygon(data = synmap_long, aes(x = x, y = y, group = alignment),
               fill = "#107AB0", alpha = .6) +
  # add prel fissions
  geom_rect(data = fissions, aes(xmin = target_min_cum, xmax = target_max_cum,
                                 ymin = ypos - 2, ymax = ypos + 2),
            color = "red") +
  geom_curve(data = fusions, aes(x = target_pos_1_cum, xend = target_pos_2_cum,
                                 y = ypos + 1, yend = ypos + 1),
             curvature = +.1) +
  geom_text(data = genomes[which(genomes$genome == "target"),], aes(x = cumulative_start + length /2,
                                                                   y = ypos + 2, label = chrom),
            angle = -60, size = 3, hjust = 1, vjust = 1) +
  geom_text(data = genomes[which(genomes$genome == "query"),], aes(x = cumulative_start + length /2,
                                                                    y = ypos - 1, label = chrom),
            angle = 60, size = 3, hjust = 1, vjust = 1) +
  theme_void() +
  ylim(c(min(genomes$ypos) - 5, max(genomes$ypos) + 5))

ggsave(paste0(outpref,"_syntenyplot_w_fissions_and_fusions.pdf"),width = 10, height = 5)
 
