p <- ggplot(pp, aes(x = V3, y = V4, color = sp)) +
  geom_point(size = 1) +  # 使用三维点图
  labs(title = "PCA through GEA loci",
       x = "PC1(63.9%)",
       y = "PC2(20.1%)") +
  theme_classic()+
  scale_color_manual(values = c("#f47c7c","#f7f48b","#70a1d7", "#a1de93","#d7a1d7","#f7c77c"))

p


# 使用 scatterplot3d 绘制三维散点图
scatterplot3d(pp$V3, pp$V4,pp$V5, 
              pch = 19,  # 点的形状
              color = pp$color,  # 根据物种着色
              main = "PCA through Polyloci",
              xlab = "PC1(38.1%)", ylab = "PC2(24.6%)", zlab = "PC3(6.16%)")

# 添加图例
legend("bottom", legend = unique(pp$sp), 
       pch = 19, col = unique(c))
