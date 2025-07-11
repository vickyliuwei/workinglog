library(ggplot2)

df <- data.frame(
  group = rep(c("GEA_loci", "Polymorphic_loci"),each = 6),
  category = rep(c("cds_variant", "5_prime_UTR_variant", "3_prime_UTR_variant","non_coding_transcript_exon_variant","intron_variant","inter_geneic_region"), 2),
  value = c(0.18,0.02,0.03,0.003,0.42,0.35,0.12,0.03,0.03,0.01,0.23,0.58)
)
# 绘制重叠条形图
pdf("snpeff.pdf",width = 10,height = 5)
ggplot(df, aes(x=group, y=value, fill=category)) +
  geom_bar(stat="identity", position="stack",width = 0.5) +
  theme_classic(base_size = 18) +
  coord_flip() +
  scale_fill_manual(values = c("#f47c7c","#f7f48b","#70a1d7", "#a1de93","#d7a1d7","#f7c77c"))
dev.off()
