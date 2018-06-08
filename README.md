# 更快的Elastic-Net|fasterElasticNet

鉴于CRAN目前的Elastic-Net算法速度较慢，因此我们考虑对其算法进行优化，运用givens变换以及前向/后向更新等方法在保证精度的同时对速度进行优化，并借助RcppArmadillo线性矩阵计算库对速度做进一步的提升，优化后速度有较大幅度的提升。

***
### Info
Version:1.0.0

***

### installation
- install `devtools`
```
install.packages(devtools)
library(devtools)
```
- install **Rcpp** and **RcppArmadillo**
```
install.packages(c('Rcpp', 'RcppArmadillo'))
```
- install **fasterElasticNet**
```
devtools::install_github('CUFESAM/Elastic-Net')
```

***

### 更新日志|Update log
#### v1.0.0
- Creat Package and upload.