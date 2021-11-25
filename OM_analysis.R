### ------------------------------------------------------------------------ ###
### analysis of OMs ####
### ------------------------------------------------------------------------ ###
library(tidyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(foreach)

### ------------------------------------------------------------------------ ###
### cod - MSY by OM ####
### ------------------------------------------------------------------------ ###
OMs_cod <- c("baseline", "rec_higher", "M_dd", "M_no_migration")
names(OMs_cod) <- OMs_cod

MSY_runs <- foreach(OM = OMs_cod, .combine = bind_rows) %do% {
  #browser()
  MSY_runs_i <- read.csv(paste0("input/cod.27.47d20/", OM, 
                                "/1000_100/MSY_trace.csv"))
  MSY_runs_i$OM <- OM
  return(MSY_runs_i)
}
MSY_runs <- MSY_runs %>%
  group_by(OM) %>%
  mutate(MSY = ifelse(catch == max(catch), TRUE, FALSE)) %>%
  mutate(catch = catch/1000, ssb = ssb/1000, rec = rec/1000) %>%
  pivot_longer(c(catch, ssb, rec)) %>%
  mutate(name = factor(name, levels = c("catch", "ssb", "rec"), 
                       labels = c("Catch [1000t]", "SSB [1000t]",
                                  "Recruitment [millions]")),
         OM = factor(OM, 
                     levels = c("baseline", "rec_higher", "M_dd", 
                                "M_no_migration"),
                     labels = c("baseline", "R: higher", "M: density dep.", 
                                "M: no migration")))
MSY_runs %>%
  ggplot(aes(x = Ftrgt, y = value)) +
  geom_vline(data = MSY_runs %>% filter(MSY == TRUE),
             aes(xintercept = Ftrgt),
             size = 0.5, colour = "blue") +
  geom_point(size = 0.3) +
  geom_line(size = 0.3) + 
  # stat_smooth(aes(alpha = "loess smoother"), size = 0.3,
  #             method = "loess", se = FALSE, span = 0.2, n = 100, 
  #             show.legend = TRUE) + 
  # scale_alpha_manual("", values = 1) +
  facet_grid(name ~ OM, scales = "free_y", switch = "y") +
  labs(x = "F (ages 2-4)") + 
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text = element_text(size = 8),
        strip.switch.pad.grid = unit(0, "pt"),
        axis.title.y = element_blank())
ggsave(filename = "output/plots/OM/OM_cod_MSY.png", 
       width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_cod_MSY.pdf",
       width = 17, height = 12, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### ple - MSY by OM ####
### ------------------------------------------------------------------------ ###
OMs_ple <- c("baseline", "M_high", "M_low", "M_Gislason",
             "no_discards_not_hidden", "rec_no_AC")

MSY_runs <- foreach(OM = OMs_ple, .combine = bind_rows) %do% {
  #browser()
  MSY_runs_i <- read.csv(paste0("input/ple.27.7e/", OM, 
                                "/1000_100/MSY_trace.csv"))
  MSY_runs_i$OM <- OM
  return(MSY_runs_i)
}
MSY_runs <- MSY_runs %>%
  group_by(OM) %>%
  mutate(MSY = ifelse(catch == max(catch), TRUE, FALSE)) %>%
  mutate(catch = catch/1000, ssb = ssb/1000, rec = rec/1000) %>%
  pivot_longer(c(catch, ssb, rec)) %>%
  mutate(name = factor(name, levels = c("catch", "ssb", "rec"), 
                       labels = c("Catch [1000t]", "SSB [1000t]",
                                  "Recruitment [millions]")),
         OM = factor(OM, 
                     levels = c("baseline", "M_high", "M_low", "M_Gislason",
                                "no_discards_not_hidden", "rec_no_AC"),
                     labels = c("baseline", "M: high", "M: low", "M: Gislason",
                                "Catch: no discards", "R: no AC")))
MSY_runs %>%
  ggplot(aes(x = Ftrgt, y = value)) +
  geom_vline(data = MSY_runs %>% filter(MSY == TRUE),
             aes(xintercept = Ftrgt),
             size = 0.5, colour = "blue") +
  geom_point(size = 0.3) +
  geom_line(size = 0.3) + 
  # stat_smooth(aes(alpha = "loess smoother"), size = 0.3,
  #             method = "loess", se = FALSE, span = 0.2, n = 100, 
  #             show.legend = TRUE) + 
  # scale_alpha_manual("", values = 1) +
  facet_grid(name ~ OM, scales = "free_y", switch = "y") +
  labs(x = "F (ages 3-6)") + 
  #scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text = element_text(size = 8),
        strip.switch.pad.grid = unit(0, "pt"),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(8, "pt"))
ggsave(filename = "output/plots/OM/OM_ple_MSY.png", 
       width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_ple_MSY.pdf",
       width = 17, height = 12, units = "cm", dpi = 600)



