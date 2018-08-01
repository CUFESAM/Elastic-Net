# 更快的Elastic-Net|fasterElasticNet

鉴于CRAN目前的Elastic-Net算法速度较慢，因此我们考虑对其算法进行优化，运用givens变换以及前向/后向更新等方法在保证精度的同时对速度进行优化，并借助RcppArmadillo线性矩阵计算库对速度做进一步的提升，优化后速度有较大幅度的提升。

由于部分Mac OS X系统中自带的clang无法支持openmp,这个分支中是去除掉makevars中对openmp的支持.

***
### Info
*Version*: 1.1.2

*Date*: 2018-07.19

*Author*: Jingyi Ma, Qiuhong Lai, Linyu Zuo, Xinyuan Yang, Xiao Liu, Yu Bai, Meng Su, Yi Yang

*Maintainer*: Jingyi Ma <jingyima@163.com>

*License*: GPL (>= 2)

*Imports*: Rcpp (>= 0.12.16)

*LinkingTo*: Rcpp, RcppArmadillo

*NeedsCompilation*: yes

***

### installation
- install `devtools`
```
install.packages(devtools)
```
- install `Rcpp` and `RcppArmadillo`
```
install.packages(c('Rcpp', 'RcppArmadillo'))
```
- install `fasterElasticNet`
```
devtools::install_github('CUFESAM/Elastic-Net')
```

- install `fasterElasticNet` without openmp supporting *Usually under clang with xcode*
```
devtools::install_github('CUFESAM/Elastic-Net', ref = 'MacOSX')
```

***

### 更新日志|Update log
#### v1.1.2|18.7.19
- Add a **predict function**.
- Update Makevars

#### v1.1.1|18.6.18
- Update Rd, edit functions' name.

#### v1.1.0
- Adding Dataset 'Housing'.

#### v1.0.0
- Creat Package and upload.
