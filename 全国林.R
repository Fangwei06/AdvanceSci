Distance decay relationship (DDR)
地理距离衰减关系
2024年-4月-18日


library(vegan)   # 加载vegan包，常用于生态学数据分析
library(geosphere)  # 加载geosphere包，用于计算地理距离
library(ggplot2)  # 加载ggplot2包，用于数据可视化
library(RColorBrewer)  # 加载RColorBrewer包，用于颜色选择

# 读取OTU数据
otu_ba <- read.csv("otutab真菌.csv", row.names=1)
# 读取分类因子数据
fa <- read.csv("metadata.csv", row.names=1)

# 按分类因子选择数据子集
se_ar <- subset(fa, cls == "BF")
se_fo <- subset(fa, cls == "CF")
se_de <- subset(fa, cls == "MF")
#se_gr <- subset(fa, cls == "Grass")
#se_we <- subset(fa, cls == "Wetland")

# 计算各子集的Bray-Curtis距离矩阵
beta_ba <- as.matrix(vegdist(otu_ba))

# 分别计算每个生态类别的Bray-Curtis距离和地理距离，并存为数据框
beta <- as.dist(beta_ba[rownames(se_ar), rownames(se_ar)])
geo.dist <- as.dist(distm(fa[rownames(se_ar), 3:4]))
df1 <- as.data.frame(cbind(bray = beta, dist = geo.dist))
df1$p <- "BF"

beta <- as.dist(beta_ba[rownames(se_fo), rownames(se_fo)])
geo.dist <- as.dist(distm(fa[rownames(se_fo), 3:4]))
df2 <- as.data.frame(cbind(bray = beta, dist = geo.dist))
df2$p <- "CF"

beta <- as.dist(beta_ba[rownames(se_de), rownames(se_de)])
geo.dist <- as.dist(distm(fa[rownames(se_de), 3:4]))
df3 <- as.data.frame(cbind(bray = beta, dist = geo.dist))
df3$p <- "MF"


# 合并所有分类的数据框
df <- rbind(df1, df2, df3)
df$p <- factor(df$p, levels = c("BF", "CF", "MF"))

# 配置颜色
col1 <- brewer.pal(8, "Set1")

# 绘制距离-相似性关系图并保存为PDF

pdf("DDRt测试颜色.pdf",width=8,height=2.5)
p<-ggplot(df,aes(dist/10000,100-bray*100,color=p))
p + 
  geom_point(size=1.5,alpha=0.9)+
  stat_smooth(method ='lm',color="black")+
  facet_wrap(~p,nrow=1)+
  scale_color_manual(values=col1)+
  xlab("Geographic distance (10 km)")+ylab("Community similarity (%)")+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()




# 加载必要的包
library(ggpmisc)  # 加载ggpmisc包，用于在图表中添加回归模型的R²和p值

# 配置颜色
col1 <- brewer.pal(8, "Set1")

# 绘制距离-相似性关系图并保存为PDF
pdf("细菌测试DDRt_with_R2_Pvalue.pdf", width=8, height=2.5)

p <- ggplot(df, aes(dist/10000, 100 - bray * 100, color = p)) +
  geom_point(size = 1.5, alpha = 0.9) +
  stat_smooth(method = 'lm', color = "black") +
  facet_wrap(~p, nrow = 1) +
  scale_color_manual(values = col1) +
  xlab("Geographic distance (10 km)") + 
  ylab("Community similarity (%)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# 添加R²和p值到图表中
p <- p + stat_poly_eq(aes(label = paste( 
                                        ..rr.label.., 
                                        ..p.value.label.., sep = "~~~")), 
                      formula = y ~ x, parse = TRUE, label.x = "right", size = 3)

print(p)

dev.off()

######加上斜率#############
# 加载必要的包
# 加载必要的包

# 计算并输出每个类别的斜率值
# 创建用于计算的数据框
df_slope <- df
df_slope$dist_10km <- df_slope$dist / 10000  # 转换为10公里单位，与图中一致
df_slope$similarity <- 100 - df_slope$bray * 100  # 计算群落相似性百分比

# 按类别分组计算斜率
slopes <- by(df_slope, df_slope$p, function(sub_df) {
  lm_fit <- lm(similarity ~ dist_10km, data = sub_df)
  slope <- coef(lm_fit)[["dist_10km"]]
  return(slope)
})

# 转换为数据框格式
slope_df <- data.frame(
  Category = names(slopes),
  Slope = as.numeric(slopes)
)

# 输出到控制台
cat("每个类别的斜率值（每增加10公里，群落相似性变化的百分比）：\n")
print(slope_df)

# 可选：保存到CSV文件
write.csv(slope_df, "DDR_真菌slope_values.csv", row.names = FALSE)







############原始代码#################


library(vegan)
library(geosphere) 
library(ggplot2)
library(RColorBrewer)

otu_ba<-read.csv("ASV_ba_rep.csv",row.names=1)
fa<-read.csv("factor_rep.csv",row.names=1)

se_ar<-subset(fa,cls=="Agriculture")
se_fo<-subset(fa,cls=="Forest")
se_de<-subset(fa,cls=="Desert")
se_gr<-subset(fa,cls=="Grass")
se_we<-subset(fa,cls=="Wetland")

beta_ba<-as.matrix(vegdist(otu_ba))
beta<-as.dist(beta_ba[rownames(se_ar),rownames(se_ar)])
geo.dist<-as.dist(distm(fa[rownames(se_ar),3:4]))
df1<-as.data.frame(cbind(bray=beta,dist=geo.dist))
df1$p<-"Agriculture"
beta<-as.dist(beta_ba[rownames(se_fo),rownames(se_fo)])
geo.dist<-as.dist(distm(fa[rownames(se_fo),3:4]))
df2<-as.data.frame(cbind(bray=beta,dist=geo.dist))
df2$p<-"Forest"
beta<-as.dist(beta_ba[rownames(se_we),rownames(se_we)])
geo.dist<-as.dist(distm(fa[rownames(se_we),3:4]))
df3<-as.data.frame(cbind(bray=beta,dist=geo.dist))
df3$p<-"Wetland"
beta<-as.dist(beta_ba[rownames(se_de),rownames(se_de)])
geo.dist<-as.dist(distm(fa[rownames(se_de),3:4]))
df4<-as.data.frame(cbind(bray=beta,dist=geo.dist))
df4$p<-"Desert"
beta<-as.dist(beta_ba[rownames(se_gr),rownames(se_gr)])
geo.dist<-as.dist(distm(fa[rownames(se_gr),3:4]))
df5<-as.data.frame(cbind(bray=beta,dist=geo.dist))
df5$p<-"Grass"
df<-rbind(df1,df2,df3,df4,df5)
df$p<-factor(df$p,levels=c("Agriculture","Forest","Wetland","Grass","Desert"))

col1<-brewer.pal(8,"Dark2")
pdf("DDRt.pdf",width=4.4,height=10)
p<-ggplot(df,aes(dist/10000,100-bray*100,color=p))
p + 
  geom_point(size=1.5,alpha=0.3)+
  stat_smooth(method ='lm',color="black")+
  facet_wrap(~p,nrow=5)+
  scale_color_manual(values=col1)+
  xlab("Geographic distance (10 km)")+ylab("Community similarity (%)")+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()

