library(ggplot2)
library(dplyr)

# Read and combine all chromosome results
files <- list.files("~/pwas", pattern="schizophrenia\\.chr[0-9]+\\.dat$", full.names=TRUE)
dat_list <- lapply(files, function(f) {
  d <- tryCatch(read.table(f, header=TRUE, as.is=TRUE, sep="\t"), error=function(e) NULL)
  if (!is.null(d)) d[!is.na(d$TWAS.P), ] else NULL
})
dat <- do.call(rbind, dat_list[!sapply(dat_list, is.null)])

# Force numeric columns (some may read as character if any row has literal "NA")
dat$TWAS.Z <- suppressWarnings(as.numeric(dat$TWAS.Z))
dat$TWAS.P <- suppressWarnings(as.numeric(dat$TWAS.P))
dat$CHR    <- suppressWarnings(as.integer(dat$CHR))
dat$P0     <- suppressWarnings(as.numeric(dat$P0))
dat$P1     <- suppressWarnings(as.numeric(dat$P1))

# Drop rows with missing TWAS result
dat <- dat[!is.na(dat$TWAS.Z) & !is.na(dat$TWAS.P), ]

# Gene name only (after the last dot)
dat$GENE <- sub(".*\\.", "", dat$ID)

# Midpoint of each feature's genomic coordinates
dat$MID <- (dat$P0 + dat$P1) / 2

# Compute cumulative x positions so chromosomes sit side by side
chr_info <- dat %>%
  group_by(CHR) %>%
  summarise(max_pos = max(MID)) %>%
  arrange(CHR)
chr_info$offset <- c(0, cumsum(head(chr_info$max_pos, -1)))
dat <- merge(dat, chr_info[, c("CHR", "offset")], by="CHR")
dat$x_pos <- dat$MID + dat$offset
chr_info$label_pos <- chr_info$offset + chr_info$max_pos / 2

# Alternating colours for chromosomes
dat$col_group <- ifelse(dat$CHR %% 2 == 0, "A", "B")

# Bonferroni Z threshold (two-sided, N = 1761 features in the ROSMAP panel)
n_feat <- 1761
bonf_z <- qnorm(0.05 / (2 * n_feat), lower.tail = FALSE)

# Significant hits for labelling
sig <- dat[abs(dat$TWAS.Z) > bonf_z, ]

p <- ggplot(dat, aes(x = x_pos, y = TWAS.Z, colour = col_group)) +
  geom_point(size = 1.2, alpha = 0.7) +
  geom_hline(yintercept = c(bonf_z, -bonf_z), colour = "red",
             linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
  geom_point(data = sig, aes(x = x_pos, y = TWAS.Z),
             colour = "firebrick", size = 2.5) +
  geom_text(data = sig, aes(x = x_pos, y = TWAS.Z, label = GENE),
            colour = "black", size = 3, vjust = -0.8, fontface = "bold") +
  scale_colour_manual(values = c("A" = "#4393C3", "B" = "#2166AC")) +
  scale_x_continuous(breaks = chr_info$label_pos, labels = chr_info$CHR,
                     expand = c(0.01, 0)) +
  labs(
    x = "Chromosome",
    y = "PWAS Z-score",
    title = "Genome-wide PWAS: Schizophrenia risk",
    subtitle = paste0("ROSMAP dorsolateral prefrontal cortex proteome (n=376)  |  ",
                      "Red dashed lines = Bonferroni threshold (p < 0.05/", n_feat, ")  |  ",
                      "Positive Z = higher protein abundance associated with greater risk")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

out <- "~/pwas/miami_plot.png"
ggsave(out, p, width = 14, height = 6, dpi = 150)
cat("Plot saved to", out, "\n")
