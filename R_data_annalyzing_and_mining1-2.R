
# 赋值，=（=赋值）、<-(向左赋值,最常用)、->(向右赋值)
# "alt" + "-"为 "<- "
# " ctrl"+"L"为清除控制台

#  查看帮助

?mean
?sd
help(sort)
help(class)

a=1
a
b <- 2
b
3->c
c

ab <- "hello,world" 
ab
a1 <- c(1,2,3,4,5,6,7,8,9)
a1
b1 <- c("asd","fgh","zxc")
b1

# 缺失值，用NA表示

x <- c(1,2,3,4,5,6,7,8,9,NA)
x

#  判断x是否含有缺失值

is.na(x)  

#  数据类型

v <- TRUE
class(v)    # class()判断数据类型

v <- 23.5
class(v)
z <- 23
class(z)
y <- "a"
class(y)

#  向量及基本运算

euro
rivers
state.abb
state.name
state.area
x <- c(1,2,3,4,5)
x
y <- c("one","two","three")
z <- c(TRUE,FALSE,T,F)
c(1:100)
seq(1,100)
seq(1,100,by=2)
seq(1,100,length.out = 10)
rep(2,5)
rep(x,10)
rep(x,each=5)
a <- c(1,2,"one")
mode(a)  # 查看数据类型
typeof(a)   # 查看数据类型
mode(z)

a=2
b=3
c="hello,world"
d=TURE
x=c(1,2,3,4,5)
y=6:10
x*2+y
x[x>3]
rep(x,c(2,4,6,1,3))
x <- c(1:100)
length(x)
# vector numberic index
x[1]
x[0]
x[-19]
x[4:18]
x[c(1,23,45,67,89)]
x[c(11,11,23,234,5,90,2)]
x[c(-2,3,4)]     #  这里是逻辑问题，不能既排除2又只显示某一个值 #
# vector logic index
y <- c(1:10)
y[c(T,F,T,T,F,F,T,T,T,F,T)]
y[c(T)]
y[c(T,F)]
y[c(T,F,T)]
y[y>5]
y[y>5 & y<9]

#charactor vector
z <- c("one","two","three","four","five")
"one" %in% z
z[z %in% c("one","two")]

# vector name index
names(y) <- c("one","two","three","four","five","six","seven","eight","nine","ten")
euro
euro["ATS"]
names(euro)

#change vector value
x <- 1:100
x
x[101] <- 101
x
x[9] <- 10
x
v <- c(1,2,3)
v
v[c(4,5,6)] <- c(4,5,6)
v
v[20]=4
v
append(v,99,after = 5)   #  在第5位数后加99
v


c <- 7:11
c
c1 <- c(c,12,13,14,15)   #  在c的末尾加上12、13、14、15
c1
c2 <- c(1,2,3,4,5,6,c)   #  在c的末尾加上1:6
c2

#   基本运算与逻辑

x <- 1:10
x+1
x-3
x <- x+1
y <- seq(1,100,length.out = 10)
x+y
x*y
x**y
y%%x     #  求余
y%/%x    #  整除
sin(x)
cos(x)
tan(x)


a <- 1:10
a
b <- 2:11
b
a>b
a+1>b
a<b
a+1==b
a!=b

a <- c(1,2,1,0,3,-1,-2)
a>1|a<0
a>1||a<0
a>1&a<0
a>1&&a<0

x <- -5:5
abs(x)
sqrt(x)
log(16,base = 2)
log10(10)
log(16)
exp(x)  #  e的n次方
ceiling (c(-2.3,3.1415))  #  向右去最接近整数
floor(c(-2.3,3.1415))     #  向左去最接近整数
trunc(c(-2.3,3.1415))     #  取整数部分
round (c(-0.618,3.1415),digits=2) #  四舍五入
signif (c(-0.0618,3.1415,10.1),digits=3)    #  有效数字  

#  向量的运算

a1 <- c(1,3,-1,-2,-5)
a2 <- c(-1,-3,-1,2,5)
a1+2
a2+2
a1*2
a2*2
a1-2
a2-2
a1/2
a2/2
a1+a2
a1*a2
a1-a2
a1/a2
a1>-1
a1==a2
a1%%2==2

x1 <- c(1,3,-1,-2,-5,6)
x2 <- c(-1,-2)
x1+x2
x1-x2



#  矩阵

iris3
state.x77
heatmap(state.x77)
#Creating Matrices
x <- 1:20
m <- matrix(x,nrow = 4,ncol = 5)
m <- matrix(1:20,4,5)
matrix(x,nrow=4,ncol=6)
matrix(x,4,4)
matrix(x,3,3)
m <- matrix(x,nrow = 4,ncol = 5,byrow = TRUE)
rnames <- c("R1","R2","R3","R4")
cnames <- c("C1","C2","C3","C4","C5")
dimnames(m)=list (rnames,cnames)


y <- matrix(1:20, nrow=5, ncol=4)
y
cells    <- c(1,26,24,68)
rnames   <- c("R1", "R2")
cnames   <- c("C1", "C2") 
mymatrix <- matrix(cells, nrow=2, ncol=2, byrow=TRUE,
                   dimnames=list(rnames, cnames)) 
mymatrix
mymatrix <- matrix(cells, nrow=2, ncol=2, byrow=FALSE,
                   dimnames=list(rnames, cnames))
mymatrix


#  another method
x <- 1:20
dim(x) <- c(4,5)
x

a <- 1:20
m <- array(a,dim = c(4,5))
m

m <- matrix(x,nrow = 4,ncol = 5)
m[1,2]
m[1,c(2,3,4)]
m[c(2,4),c(2,3)]
m[2,]
m[,2]
m[2]
m[-1,2]
m[-1,]


#Using matrix names

x <- 1:20
m <- matrix(x,nrow = 4,ncol = 5,byrow = TRUE)
rnames <- c("R1","R2","R3","R4")
cnames <- c("C1","C2","C3","C4","C5")
dimnames(m)=list (rnames,cnames)
m
m["R1","C2"]
m["R1",]
m[,c("C1","C2")]

state.x77
state.x77[,"Income"]
state.x77["Alabama",]

#  向量计算

x <- 1:20
m <- matrix(x,nrow = 4,ncol = 5,byrow = TRUE)

m+1
m*2
m+m
n <- matrix(1:20,5,4)
m+n
colSums(m)
rowSums(m)
colMeans(m)
rowMeans(m)

n <- matrix (1:9,3,3)
t <- matrix (2:10,3,3)
n
t
n*t
n%*%t
diag(n)    #  对角线
diag(t)


#   数组

dim1 <- c("A1", "A2")
dim2 <- c("B1", "B2", "B3")
dim3 <- c("C1", "C2", "C3", "C4")
z <- array(1:24, c(2,3,4), dimnames=list(dim1, dim2, dim3))
z

x <- 1:20
dim(x) <- c(2,2,5)
dim1 <- c("A1", "A2")
dim2 <- c("B1", "B2", "B3")
dim3 <- c("C1", "C2", "C3", "C4")
z <- array(1:24, c(2,3,4), dimnames=list(dim1, dim2, dim3))
z


Titanic


#  数据框

patientID <- c(1, 2, 3, 4)
age <- c(25, 34, 28, 52)
diabetes <- c("Type1", "Type2", "Type1", "Type1")
status <- c("Poor", "Improved", "Excellent", "Poor")
patientdata <- data.frame(patientID, age, diabetes, status)
patientdata
colnames(patientdata) <- c("ID","Age1")
patientID


iris
mtcars
rock
state <- data.frame(state.name,state.abb,state.region,state.x77)
state[1]
state[c(2,4)]
state[,"state.abb"]
state["Alabama",]

patientdata[1:2]
patientdata[c("diabetes","status")]
patientdata$age
patientID <- c(1, 2, 3, 4)
age <- c(25, 34, 28, 52)
diabetes <- c("Type1", "Type2", "Type1", "Type1")
status <- c("Poor", "Improved", "Excellent", "Poor")
patientdata <- data.frame(patientID, age, diabetes, status)
patientdata

women
women$height
plot(women$height,women$weight)
lm (weights~height,data = women)
attach(mtcars)
mpg
hp
detach(mtcars)
with(mtcars,{mpg})




#  列表

state.center
a <- 1:20
b <- matrix(1:24,4,6)
c=mtcars
d <- "This is a test list"
mlist <- list(a,b,c,d)
mlist <- list(first=a,second=b,third=c,fourth=d)

#List
mlist[1]
mlist[c(1,4)]
state.center[c("x","y")]
mlist$first
state.center$x
mlist[[1]]
class(mlist[2])
class(mlist[[2]])
mlist[[5]] <- iris
mlist
mlist[5] <- NULL
mlist[[5]] <- NULL


g <- "My First List"
h <- c(25, 26, 18, 39)
j <- matrix(1:10, nrow=5)
k <- c("one", "two", "three")
mylist <- list(title=g, ages=h, j, k)
mylist


#   因子

patientID <- c(1, 2, 3, 4)
age <- c(25, 34, 28, 52)
diabetes <- c("Type1", "Type2", "Type1", "Type1")
status <- c("Poor", "Improved", "Excellent", "Poor")
diabetes <- factor(diabetes)
status <- factor(status, order=TRUE)
patientdata <- data.frame(patientID, age, diabetes, status)
str(patientdata)                               
summary(patientdata)

mtcars
mtcars$cyl
table(mtcars$cyl)
f <- factor(c("red","red","green","red","blue","green","blue","blue"))
week <- factor(c("Mon","Fri","Thu","Wed","Mon","Fri","Sun"))
week <- factor(c("Mon","Fri","Thu","Wed","Mon","Fri","Sun"),order = TRUE,
               levels = c("Mon","Tue","Wed","Thu","Fri","Sat","Sun"))
fcyl <- factor(mtcars$cyl)
plot(mtcars$cyl)
plot(fcyl)
num <- c(1:100)
cut (num,c(seq(0,100,10)))
state.division
state.region
is.factor(state.division)  #  判断是否为因子


#   读写数据

x <- read.table ("input.csv",sep=",",header = T)
head(x)   #  只看x前6行，常用
tail(x)   #  只看x后6行
head(x,10)   #  只看x前10行
x <- read.table ("input 1.txt",sep=",",header = T,skip = 5)    #  跳5行
x <- read.table ("input.csv",sep=",",header = T,nrows = 100)    #  只看100行
x <- read.table ("input.csv",sep=",",header = T,skip = 5,nrows = 100)

x <- read.csv ("input.csv",header = T)
x <- read.csv ("input.csv",sep = ",",header = T)

write.table (x,file=newfile.txt)
write.table (x,file=newfile.csv,sep="\t")
write.table (x,file=newfile.csv,sep="\t",quote=FALSE,append=FALSE,na="NA")
write.table (x,file=gzfile (newfile.csv.gz),sep="\t",
             quote=FALSE,append=FALSE,na="NA")


install.packages("xlsx")   #  安装R包
library(xlsx)     #  加载R包
rdata <- read.xlsx("data.xlsx",sheetIndex = 1,startRow = 1,endRow = 100)
write.xlsx(rdata,file = "rdata.xlsx",sheetName = "Sheet 1",append = F)
help(package="xlsx")

#  排序

mtcars
sractm <- t(mtcars)   #  t将表格方向改变
letters
rev(letters)   #  rev排列倒换 
women[rev(rownames(women)),]
transform(women, height = height*2.54)
transform(women, cm = height*2.54)
#sort order rank
sort(rivers)   #  sort从小到大排列
sort(state.name)
rev(sort(rivers))
mtcars[sort(rownames(mtcars)),]
sort(rivers)
order (rivers)  #  排列的排名
mtcars[order(mtcars$mpg),]
mtcars[order(-mtcars$mpg),]    #   -表示反向
order(mtcars$mpg,mtcars$disp)

#  常用统计分析

attach(mtcars)
head(mtcars)
mean(mpg)
median(mpg)
quantile(mpg)    #  分位数
quantile(mpg,probs = c(0.25,0.75))
IQR(mpg)    #  四分位差

#  众数
#  1个众数

table(gear)
which.max(table(gear))

#  多个众数

x <- table(carb)
tmp <- x
tmp.max <- max(tmp)
which(tmp==tmp.max)

range(mpg)
diff(range(mpg))
sd(mpg)
var(sd)
cv <- sd/mean
scale(mpg)
detach(mtcars)

#  如果不用attach,每一次都要说明操作对象

mtcars
mean(mtcars$mpg)
var(mtcars$mpg)

WorldPhones
worldphones <- as.data.frame(WorldPhones)
rs <- rowSums(worldphones)
cm <- colMeans(worldphones)
total <- cbind(worldphones,Total=rs)
rbind(total,cm)

#apply series

apply(WorldPhones,MARGIN = 1,FUN = sum)
apply(WorldPhones,MARGIN = 2,FUN = mean)

#lapply and sapply
state.center
lapply(state.center,FUN = length)
sapply(state.center,FUN = length)
#tapply
state.name
state.division
tapply(state.name,state.division,length)

#scale and center
state.x77
heatmap(state.x77)
x <- scale(state.x77)
heatmap(x)


x <- c(1,2,3,6,3)
mean(x)
x - mean(x)
sd(x)
(x -mean(x))/sd(x)
scale(x)



#   偏度和峰度


#  自己写函数

mystats <- function(x,na.omit=FALSE) {
  if(na.omit) 
    x <- x[!is.na(x)]
  m <- mean(x)
  n <- length(x)
  s <- sd(x)
  skew <- sum(((x-m)/s)^3)/n     #  偏度
  kurt <- sum((((x-m)/s)^4)/n)-3    #  风度
  return(c(n=n,mean=m,stdev=s,skew=skew,kurtosis=kurt))
}

x <- 1:10
mystats(x)

#  别人的包


install.packages("agricolae")    #  安装包
library(agricolae)    #加载包
skewness(x)
kurtosis(x)

#   按照公式，我觉得自己写保险




























































