# 更快的Elastic-Net|fasterElasticNet

鉴于CRAN目前的Elastic-Net算法速度较慢，因此我们考虑对其算法进行优化，运用givens变换以及前向/后向更新等方法在保证精度的同时对速度进行优化，并借助RcppArmadillo线性矩阵计算库对速度做进一步的提升，优化后速度有较大幅度的提升。

***
### Info
*Version*: 1.0.0

*Date*: 2018-05-30

*Author*: Jingyi Ma, Qiuhong Lai, Linyu Zuo, Xinyuan Yang, Xiao Liu, Yu Bai, Meng Su, Yi Yang

*Maintainer*: Jingyi Ma <jingyima@163.com>

　　　　&nbsp;　Linyu Zuo <zuozhe5959@gmail.com>

*License*: GPL (>= 2)

*Imports*: Rcpp (>= 0.12.16)

*LinkingTo*: Rcpp, RcppArmadillo

*NeedsCompilation*: yes

***

### installation
- install `devtools`
```
install.packages(devtools)
library(devtools)
```
- install `Rcpp` and `RcppArmadillo`
```
install.packages(c('Rcpp', 'RcppArmadillo'))
```
- install `fasterElasticNet`
```
devtools::install_github('CUFESAM/Elastic-Net')
```

***

### 更新日志|Update log
#### v1.0.0
- Creat Package and upload.