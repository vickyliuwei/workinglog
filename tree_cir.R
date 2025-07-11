library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
max_length <- 0.2  # 设定最大显示长度
tree_capped <- tree
tree_capped$edge.length[tree_capped$edge.length > max_length] <- max_length
n_tips <- length(tree$tip.label)     # 叶节点数量 (147)
n_nodes <- tree$Nnode                 # 内部节点数量 (144)
internal_nodes <- (n_tips + 1):(n_tips + n_nodes)
internal_data <- data.frame(
  label = tree$node.label,  # 内部节点无默认 label，需手动命名
  species = "Unknown"                       # 默认分类
)
tree <- read.tree("E:\\Pinus SNP paper\\tree\\109sample.treefile")
add <- gea_p[n,c(25,23)]
colnames(internal_data) <- colnames(add)
full_data <- rbind(add,internal_data)
p <- ggtree(tree_capped, layout = "circular") %<+% full_data
p <- p + geom_tiplab(aes(color = sp), size = 2, hjust = -0.1)
p <- p + geom_tippoint(aes(color = sp), size = 1)

p <- p + scale_color_manual(values = c("#f47c7c","#f7f48b","#70a1d7", "#a1de93","#d7a1d7","#f7c77c")) 

# 显示图形
print(p)
