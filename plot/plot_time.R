library(ggplot2)
library(tidyr)
df=data.frame(size=c("j = 2","j = 3","j = 4","j = 5"),
              gd=c(4.245, 8.496, 15.295, 31.4814),
              vN=c(1.09, 2.546,4.3030, 7.996),
              qN = c(3.307, 7.532, 16.3001, 32.790))
              


df$size <- factor(df$size, levels = df$size)

df <- gather(df, key=Optimizer, value="Time", gd, vN, qN)
ggplot(df, aes(x=size, y=Time, group=Optimizer,color=Optimizer)) +
  geom_line(size = 1, alpha=0.5)+
  geom_point()+
  xlab("Number of experts") +
  ggtitle("Time Complexity for 50 Iterations") + theme_minimal() +
  scale_color_discrete(name  ="Optimizers") 
