
# customised plotting theme
theme_meta <- 
  function(base_family = "") {
    theme(panel.background=element_rect(fill="white", colour="black", linetype="solid"),
          panel.grid = element_blank(),
          axis.text = element_text(colour="black",size=10),
          legend.key = element_rect(fill = NA))
  }

### END