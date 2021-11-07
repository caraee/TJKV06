#Figure 1: pBW / pBG

#Figure 2: p_GSIS_BG / p_GSIS_CP

#Figure 3: p_ITT_BG / p_ITT_CP

#Figure 4: p_ArgTT_BG / p_ArgTT_Ins / p_ArgTT_Gcg / p_ArgTT_GLP1

#Figure 2: (pBW | pBG | pHbA1c) / (p_GSIS_BG | p_GSIS_ratCP) / p_GSIS_CP
library("patchwork")

#need to do fasted/random fed comparison; is there a mouse C-peptide panel I can include?
#then finish checking TJKV03 models, TJKV07 ITT CP with fixed divergence
#then check BIRKO model
(Figure1<- (p_weeklyBW / p_weeklyBG) + plot_layout(guides = 'collect')+
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(family = "Arial",color="black",size=14)))

ggsave('./Figures/Figure1.svg', Figure1, width=8, height=5, units = "in", dpi=500)
ggsave('./Figures/Figure1.png', Figure1, width=8, height=5, units = "in", dpi=500)

(Figure2<- p_GSIS_BG / p_GSIS_CP  + plot_layout(guides = 'collect')+
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(family = "Arial",color="black",size=14)))

ggsave('./Figures/Figure2.svg', Figure2, width=12, height=8, units = "in", dpi=500)
ggsave('./Figures/Figure2.png', Figure2, width=12, height=8, units = "in", dpi=500)


layout <- '
AAAAAA
BBBBB#
'
(Figure3<- wrap_plots(A = p_ITT_BG, B = p_ITT_CP, design = layout) +
  plot_layout(guides = 'keep')+
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(family = "Arial",color="black",size=14)))

ggsave('./Figures/Figure3.svg', Figure3, width=8, height=8, units = "in", dpi=500)
ggsave('./Figures/Figure3.png', Figure3, width=8, height=8, units = "in", dpi=500)



layout4 <- '
AAAAAAAA
BBBBBBBB
CCCCC###
'

patch<-p_ArgTT_Ins | p_ArgTT_Gcg
(Figure4<-wrap_plots(A = p_ArgTT_BG, B = patch, C = p_ArgTT_GLP1,
                     design = layout4,guides="keep") +
    # plot_layout(guides = 'keep')+
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(family = "Arial",color="black",size=14)))

ggsave('./Figures/Figure4.svg', Figure4, width=8, height=12, units = "in", dpi=500)
ggsave('./Figures/Figure4.png', Figure4, width=8, height=12, units = "in", dpi=500)

inset_element(
  p,
  left,
  bottom,
  right,
  top,
  align_to = "panel",
  on_top = TRUE,
  clip = TRUE,
  ignore_tag = FALSE
)

