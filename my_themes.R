mytheme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.y = element_text(size=15), axis.title=element_text(size=15),legend.title = element_text(size=17), legend.text = element_text(size=15), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black",fill="white"),strip.text = element_text(size=14,  face= "bold"))

my_theme1 <- theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

my_theme3 <- theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position="bottom",
                   legend.text=element_text(size=13), legend.title=element_text(size=15), legend.direction = "vertical")

my_theme4 <- theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_rect(colour="black",fill="white"),
                   axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position="right",
                   legend.text=element_text(size=13), legend.title=element_text(size=15), legend.direction = "vertical",strip.text = element_text(size=14,  face= "bold"))


