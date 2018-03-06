#this script calculates sma carrier probabilties, then plots the 95% confidence intervals and carrier probabilites
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(docopt)

'usage: calculate_carrier_probability.R <input_file> <output_directory>' -> doc
opts <- docopt(doc)
input <- opts$input_file # the sma_sample_stat.txt file
output <- opts$output_directory

sample_data<-read.delim(input)

#Jeffreys noninformative prior
sample_data<-mutate(sample_data,carrier_prob=pbeta(0.38, .5+di, .5+ri-di))

sample_data <- sample_data %>%
  mutate(low  = qbeta(.025, .5+di, .5+ri-di),
         high = qbeta(.975, .5+di,.5+ri-di),
         Carrier_Status=ifelse(low<=.38 & high>=.38,"possible",ifelse(high<.38,"likely","unlikely")))

# plot confidence intervals
ci <- gather(sample_data,"value","point",low,high)
ci_plot <- ggplot(ci,aes(x=point,y=sample,color=Carrier_Status))+
  geom_point(size=.75, alpha=.75)+geom_line()+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=10),
        axis.ticks.y = element_blank(),
        axis.line.x =element_line(linetype="solid"),
        axis.line.y =element_line(linetype="solid"),
        axis.title=element_text(size=14,face="bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10))+
  geom_vline(aes(xintercept = 0.38),linetype="dashed",color="red")+
  scale_color_manual(values=c("red","yellow3","gray60"), name = "Carrier Status")+
  labs(x="Posterior 95% Credible Interval for\n Proportion of SMN1 Reads",
       y="Sample")
ci_plot 

#plot carrier probabilities 
set.seed(3)
carrier_prob__plot<-ggplot(sample_data,aes(x=pi,y=carrier_prob,colour=Carrier_Status))+
  geom_point(size=.75, position=position_jitter(w=.020,  h=.020))+
  theme_bw()+
  theme(legend.position=c(.85,0.75),
        legend.background = element_rect(colour="black",fill="white"),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10),
        axis.text = element_text(size=10),
        axis.text.x = element_text(colour = "black",angle=90),
        axis.title=element_text(size=16,face="bold"))+
  geom_vline(aes(xintercept = 0.38),linetype="dashed",color="red")+
  scale_y_continuous(expand=c(0,0),limits=c(-0.025,1.09),breaks=pretty_breaks(n=10))+
  scale_x_continuous(expand=c(0,0),limits=c(-0.025,1.09),breaks=pretty_breaks(n=10))+
  scale_color_manual(values=c("red","yellow3","gray45"))+
  labs(x="Observed Proportion of SMN1 Reads",y="SMA Carrier Probability", color="Carrier Status")
set.seed(3)
carrier_prob__plot

ggsave(ci_plot,filename=paste(output,"confidence_intervals_plot.tiff",sep="/"),dpi=300, width=7.5,height=3.5,units="in")
ggsave(carrier_prob__plot,filename=paste(output,"carrier_probabilities_plot.tiff",sep="/"),dpi=300, width=7.5,height=3.5,units="in")
