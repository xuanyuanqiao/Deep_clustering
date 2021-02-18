# Deep_clustering
### Algorithm description 
This method is based on DESC, which is an unsupervised deep learning algorithm for clustering scRNA-seq data.
The algorithm constructs a non-linear mapping function from the original scRNA-seq data space to a low-dimensional feature space by iteratively learning cluster-specific gene expression representation and cluster assignment based on a deep neural network.
Please see the published paper for detail information:https://www.nature.com/articles/s41467-020-15851-3

![Image of model](https://github.com/xuanyuanqiao/Deep_clustering/blob/main/%E7%AE%97%E6%B3%95%E6%A8%A1%E5%9E%8B%E5%9B%BE.jpg)

Requirement(version of this runing code):
    
- desc ---  2.0.3
- python --- 3.6.12
- keras --- 2.1.0
- tensorflow --- 1.7.0
- scanpy --- 1.6.0
- graphviz --- 2.40.1

## References
<a id="1">[1]</a> 
MacParland, S.A., Liu, J.C., Ma, XZ. et al. Single cell RNA sequencing of human liver reveals distinct intrahepatic macrophage populations. Nat Commun 9, 4383 (2018). https://doi.org/10.1038/s41467-018-06318-7
<a id="2">[2]</a> 
Xiangjie Li, Yafei Lyu, Jihwan Park, Jingxiao Zhang, Dwight Stambolian, Katalin Susztak, Gang Hu, Mingyao Li. Deep learning enables accurate clustering and batch effect removal in single-cell RNA-seq analysis. 2019. bioRxiv 530378; doi: https://doi.org/10.1101/530378
