library(ggplot2)
library(dplyr)
library(tidyr)

sample_data<-read.delim("ice_sma_sample_stat.txt")

#Jeffreys noninformative prior
sample_data<-mutate(sample_data,carrier_prob=pbeta(0.38, .5+di, .5+ri-di))

sample_data <- sample_data %>%
  mutate(low  = qbeta(.025, .5+di, .5+ri-di),
         high = qbeta(.975, .5+di,.5+ri-di),
         Carrier_Status=ifelse(low<=.38 & high>=.38,"possible",ifelse(high<.38,"likely","unlikely")))


ci <- gather(sample_data,"value","point",low,high)
ci_plot <- ggplot(ci,aes(x=point,y=sample,color=Carrier_Status))+
  geom_point(alpha=.75,size=3)+geom_line()+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=18),
        axis.ticks.y = element_blank(),
        axis.line.x =element_line(linetype="solid"),
        axis.line.y =element_line(linetype="solid"),
        axis.title=element_text(size=22,face="bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18),
        legend.key=element_rect(fill=NA))+
  geom_vline(aes(xintercept = 0.38),linetype="dashed",color="red")+
  scale_color_manual(values=c("red","yellow3","gray60"), name = "Carrier Status")+
  labs(x="Posterior 95% Credible Interval for\n Proportion of SMN1 Reads",
       y="Sample")
ci_plot 

