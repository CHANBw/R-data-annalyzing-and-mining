
#####################################
##              第二节             ##
#####################################
data()
class(rivers)
plot(rivers)
hist(rivers)
class(state.abb)
pie(table(state.abb))
color <- c("red","green","blue","yellow")    #  设置颜色

pie(table(state.abb),col = color)
labs <- c("NC","S","N","W")    #  设置标签
pie(table(state.abb),labels = labs)
class(state.division)
plot(mtcars$cyl,mtcars$disp)   #  点状图

plot(as.factor(mtcars$disp),mtcars$cyl)    #  箱型图,as.factor强行变因子

heatmap(state.x77)
heatmap(as.data.frame(state.x77)) 

#  该行代码无法执行，因为热图必须为矩阵，as.data.frame强行变为数据框。

fit <- lm(weight~height,data=women)
#  lm为线性模型，~表示关联
class(fit)
plot(fit)



#####################################
##              第三节             ##
#####################################
example("barplot")   #  展示例图
example("heatmap")
library(ggplot2)
example("qplot")
help("barplot")
??barplot


#####################################
##          第四节 barplot         ##
#####################################
getwd()

#  显示当前工作路径

setwd(dir = "c:/Users/wangtong/Desktop/")

#  设置工作路径，一般与你需要执行的文件处于同一路径

m <- read.csv("homo_length.csv",header = T)
m
class(m)
head(m)
x <- m[1:24,]

#  保证只有组装到染色体上的才能用来做图

x

class(x$length)
barplot(height = x$length)
barplot(height = x$length,names.arg = x$chr)

#  按染色体顺序命名

color <- sample(colours(),24,replace = F)

#  随机抽取颜色并设置，sample表抽样，colours颜色函数，
#抽24种，replace表是否放回即同一样本是否可以重复抽到

barplot(height = x$length,names.arg = x$chr,col = color)
barplot(height = x$length,names.arg = x$chr,col = rainbow(7))
barplot(height = x$length,names.arg = x$chr,col = rainbow(7),border = F)
barplot(height = x$length,names.arg = x$chr,col = rainbow(7),border = F,
        main = "Human chromosome length distribution",xlab = "Chromosome Name",
        ylab = "Chromosome Length")

#  border表示边界，main表示标题，x/ylab即x/y轴题目


#####################################
##        第5节 分组条形图         ##
#####################################
 
x <- read.csv("sv_distrubution.csv",header = T,row.names = 1)
x
barplot(x)
barplot(as.matrix(x))

#  转变为矩阵

barplot(t(as.matrix(x)))

#  转变矩阵方向

barplot(t(as.matrix(x)),col = rainbow(4))
barplot(t(as.matrix(x)),col = rainbow(4),beside = T)

#  有堆叠复试条形图变为并列复试条形图

barplot(t(as.matrix(x)),col = rainbow(4),legend.text = colnames(x))

#  加x轴标签

barplot(t(as.matrix(x)),col = rainbow(4),legend.text = colnames(x),ylim = c(0,800))

#  改变Y轴轴距

barplot(t(as.matrix(x)),col = rainbow(4),legend.text = colnames(x),ylim = c(0,800),
        main = "SV Distribution",xlab="Chromosome Number",ylab="SV Numbers")


#####################################
##        第6节 直方图             ##
#####################################

x <- read.table("H37Rv.gff",sep = "\t",header = F,skip = 7,quote = "")

# gff是一种基因注释文件，包括每个基因的大致位置，
#  CDS序列长度和位置，外显子长度、数量和位置等等，用于统计,quote=""表示里面的字符串不用“”括起来。

x <- x[x$V3=="gene",]

#  只提取第三列为gene的行，用逻辑值做索引

x <- abs(x$V5-x$V4+1)

#  计算基因长短

length(x)

#  统计又多少个基因

range(x)

#  最长基因和最短基因大小

hist(x)
hist(x,breaks = 80)

#  设置条目，即有多少个方块

hist(x,breaks = c(0,500,1000,1500,2000,2500,15000))
hist(x,breaks = 80,freq = F)
#    用概率密度画图,默认为频率
hist(x,breaks = 80,density = T)
#     是否绘制斜线
hist(rivers,density = T,breaks = 10)
?hist


#    绘制概率密度
h=hist(x,nclass=80,col="pink",xlab="Gene Length (bp)",main="Histogram of Gene Length");

#    将变量保存至h，并添上标题

h
h$counts
rug(x)

#    看分布频率

xfit<-seq(min(x),max(x),length=100)

#    将基因长度平均分为100份

yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))

#     生呈正态分布

yfit <- yfit*diff(h$mids[1:2])*length(x)

#    计算相应概率值

lines(xfit, yfit, col="blue", lwd=2)

#     绘制密度曲线，lwd表线宽。


#####################################
##        第7节 散点图             ##
#####################################
m <- read.table("prok_representative.txt",sep="\t")
head(m)
genome_size <- m[,2]
gene_number <- m[,4]
plot(genome_size,gene_number,pch=16,cex=0.4,xlab="Genome Size",ylab="Genes")

#    pch为点的形状类型，有0~25；cex为点的大小

fit <- lm(gene_number ~ genome_size)

#    建立线性模型

abline( fit,col="blue",lwd=1.8 )

#    绘制回归线

summary(fit)

#    查看统计量

summary(fit)$adj.r.squared
rr <- round( summary(fit)$adj.r.squared,2)

#    查看R的平方

summary(fit)$coefficients[1]
intercept <- round( summary(fit)$coefficients[1],2)
intercept

#    回归方程截距，round保留2为小数

summary(fit)$coefficients[2]
slope <- round( summary(fit)$coefficients[2],2)
slope

#    查看回归方程斜率

eq <- bquote( atop( "y = " * .(slope) * " x + " * .(intercept), R^2 == .(rr) ) )

#    bquote将写入内容一模一样输出或赋值，atop表输入数学表达式

text(12,6e3,eq)

#    将回归方程和R^2加入到图中，12和6e3表示方程写到哪个位置，是坐标值(12,6e3)。



#####################################
##        第8节 饼图               ##
#####################################
x <- read.csv("homo_length.csv",header = T)
x <- x[1:24,]
barplot(height = x$length,names.arg = x$chr)
x$length/sum(x$length)

#  计算每条染色体长度占总长的比例

pie(x$length/sum(x$length))

#    正式绘制饼图

m <- read.table("Species.txt")
m
x <- m[,3]
pie(x)
pie(x,col=rainbow(length(x)))

#    rainbow需要指明参数，即有多少元素需要添加颜色

paste(m[,1],m[,2],"\n",m[,3],"%" )

#    利用paste函数将每个物种的三行信息粘到一块，“\n”表换行

lbls <- paste(m[,1],m[,2],"\n",m[,3],"%" )
pie(x,col=rainbow(length(x)),labels = lbls)
pie(x,col=rainbow(length(x)),labels = lbls,radius = 1)

#    radius表示半径

pie(x,col=rainbow(length(x)),labels = lbls,radius = 1,cex=0.5)

#    cex设置文字的大小

#3D饼图

install.packages("plotrix")
library(plotrix)
pie3D(x,col=rainbow(length(x)),labels = lbls)
pie3D(x,col=rainbow(length(x)),labels = lbls,explode = 0.1)

#    explode,炸裂开，各元素分开的距离

pieplot <- pie3D(x,col=rainbow(length(x)),radius = 1,explode = 0.1)

#    将图像信息赋值给一个对象，方便处理标签

pie3D.labels(pieplot,labels = lbls,labelcex = 0.8,height = 0.1,labelrad = 1.7)

#    labelcex表示文字大小，height表示标签的高度，labelrad表标签距离饼的距离


#扇形图
fan.plot(x,col=rainbow(length(x)),labels = lbls,cex=0.8,radius = 1)



#####################################
##        第9节 箱线图             ##
#####################################
boxplot(mpg ~cyl,data = mtcars)

m <- read.csv("gene_expression.csv",header = T,row.names = 1)
head(m)
boxplot(m)

#  因为离群值太多，所以画不了

boxplot(m,outline=F)

#   去掉离群值

#   下面用公式来画

install.packages("reshape2")
library(reshape2)
x <- melt(m)

#  将数据分为两列，变量和数值

head(x)
boxplot(value ~ variable,data = x)
boxplot(value ~ variable,data = x,outline=F)

#   下面美化步骤

head(ToothGrowth)
boxplot(len ~ dose, data = ToothGrowth)

#  按剂量分组

boxplot(len ~ dose:supp, data = ToothGrowth)

#  按剂量和处理方式分组

boxplot(len ~ dose:supp, data = ToothGrowth,col = c("orange", "red"))

#  添加颜色

boxplot(len ~ dose:supp, data = ToothGrowth,col = c("orange", "red"),notch=T)

#  在中位数处开缺口

boxplot(len ~ dose:supp, data = ToothGrowth,col = c("orange", "red"),
        width=c(1,0.5,1,0.5,1,0.5))

#  每个箱体设置宽度

boxplot(len ~ dose:supp, data = ToothGrowth,col = c("orange", "red")
        ,varwidth=T)

#  箱体宽度跟随数值改变

boxplot(len ~ dose:supp, data = ToothGrowth,col = c("orange", "red")
        ,boxwex=0.5)

#  设置所有箱体宽带

boxplot(len ~ dose:supp, data = ToothGrowth,col = c("orange", "red")
        ,staplewex=0.5)

#  设置最大值和最小值（最上面和最下面两条线的宽度）

boxplot(len ~ dose:supp, data = ToothGrowth,col = c("orange", "red")
        ,sep = ":",lex.order = T)

#  sep变量名称连接方式，lex.order改变变量排列顺序


#####################################
##   第10节 par函数  设置参数      ##
#####################################

#  主要讲解绘图参数

opar <- par(no.readonly = T)
?par
par()
x <- par()
x$pch
par(pch=16)
opar <- par(no.readonly = T)
opar
par(opar)

# 一次输一节，并去修改，看函数功能

plot(women$height,women$weight)
plot(women$height,women$weight,type = "b",pch=16,col="red",lty=2,lwd=2,
     main = "Main Title",sub = "Sub Title",xlab = "Heigth",
     ylab = "Weight",xlim = c(50,100),ylim = c(100,200),cex=1.5,las=3,
     adj=0.3,bty="c",fg="blue",tck=-0.03,col.main="green",cex.lab=2)



#####################################
##        第11节 颜色              ##
#####################################
colors()
sample(colours(),24)

rgb(0, 1, 0)
hsv(0,1,0)
hcl(0,1,0)

#    以上三个函数为颜色转换，可以结合配图网站使用

rainbow(7)

#   不止7种颜色

pie(1:7,col=rainbow(7))
heat.colors(7)
pie(1:7,col=heat.colors(7))

#    暖色调

terrain.colors(7)
pie(1:7,col=terrain.colors(7))

#   素雅且土

topo.colors(7)
pie(1:7,col=topo.colors(7))

#    素雅

cm.colors(7)
pie(1:7,col=cm.colors(7))

#    由蓝到紫

gray.colors(7)
pie(1:7,col=gray.colors(7))

#    灰色


#   下面为调色板

install.packages("RColorBrewer")
library(RColorBrewer)

display.brewer.all()
#    显示调色板中多种颜色
display.brewer.pal(name = "Set1",n=9)

#    显示某一具体子集的颜色

brewer.pal.info

#    调色板信息，包括种类，是否为连续型以及红绿色盲是否可分辨

brewer.pal(name = "Set1",n=9)

#    取颜色，set1只有9种

#   下面时举例应用

color=brewer.pal(name = "Set1",n=7)
pie(1:7,col = color)



#####################################
##        第12节 热图1             ##
#####################################
#1 heatmap()
#2 gplots 包 heatmap.2()
#3 pheatmap包 pheatmap
#4 ComplexHeatmap包

#   这一节主要使用内置函数heatmap()


example("heatmap")
m <- read.csv("heatmap.csv",header = T,row.names = 1)
class(m)

#    查看数据类型

x <- as.matrix(m)

#    绘制热图数据类型需要为矩阵，强制转化为矩阵

heatmap(x)
heatmap(t(x))

#    将样品名和基因名转置

heatmap(x,col=c("green","red"))

colorRampPalette(c("red","black","green"))(nrow(x))

#    根据数据有多少行生成相应的渐变色

color <- colorRampPalette(c("red","black","green"))(nrow(x))
heatmap(x,col=color)
heatmap(x,col=color,RowSideColors = color)

#     给相应的行生成颜色跳

color1=colorRampPalette(c("red","black","green"))(ncol(x))

#    根据数据有多少列生成相应的渐变色

heatmap(x,col=color,ColSideColors = color1)

#     给相应的行生成颜色跳


#     下面是聚类相关，默认聚类

heatmap(x,col=color,Rowv = NA)
heatmap(x,col=color,Rowv = NA,Colv = NA)

#    NA表示取消相应的聚类

heatmap(state.x77)
heatmap(state.x77,scale = "col")

#####################################
##        第13节 热图2             ##
#####################################

#install.packages("gplots")
library(gplots)
m <- read.csv("heatmap.csv",header = T,row.names = 1)
x <- as.matrix(m)
heatmap.2(x)
#   用heatmap.2和上一个例子中的数据绘制热图

example("heatmap.2")

heatmap.2(x)
heatmap.2(x,key = F)
#    有无图例

heatmap.2(x,symkey = F)
heatmap.2(x,symkey = T)
#    图例颜色是否对称和连续

heatmap.2(x,symkey = T,density.info = "none")
#    图例信息

heatmap.2(x,symkey = T,trace = "none")
#    去掉热图中的竖线

heatmap.2(x,symkey = T,tracecol = "black")
#    设置竖线颜色

heatmap.2(x,dendrogram = "none")
#    是否显示聚类信息

heatmap.2(x,dendrogram = "row")
heatmap.2(x,dendrogram = "col")
#    只显示行或列的聚类信息



#    当下流行pheatmap

#install.packages("pheatmap")

library(pheatmap)
example("pheatmap")
m <- read.csv("heatmap.csv",header = T,row.names = 1)
x <- as.matrix(m)
pheatmap(x)

#    下面设置注释

colnames(x)
annotation_col <- data.frame(CellType=factor(rep(c("N1","T1","N2","T2"),each=5)))
#    样品中N和T各有两种类型，每种五个。并根据此进行分类

#    接下来就是将分类信息与行号映射起来

annotation_col
rownames(annotation_col)
#    查看行号

rownames(annotation_col) <- colnames(x)
#    完成映射

pheatmap(x,annotation_col = annotation_col)
#    显示映射结果

pheatmap(x,annotation_col = annotation_col,display_numbers = T)
#    显示每个各自中的数据

pheatmap(x,annotation_col = annotation_col,display_numbers = T,
         number_format = "%.2f")
#    保留两位有效数字

pheatmap(x,annotation_col = annotation_col,display_numbers = T,
         number_format = "%.1f",number_color = "black")
#    设置格子中字体颜色



#####################################
##        第14节 韦恩图            ##
#####################################

#install.packages("VennDiagram")


listA <- read.csv("genes_list_A.txt",header=FALSE)
A <- listA$V1
listB <- read.csv("genes_list_B.txt",header=FALSE)
B <- listB$V1
listC <- read.csv("genes_list_C.txt",header=FALSE)
C <- listC$V1
listD <- read.csv("genes_list_D.txt",header=FALSE)
D <- listD$V1
listE <- read.csv("genes_list_E.txt",header=FALSE)
E <- listE$V1
length(A);length(B);length(C);length(D);length(E)
#    查看样品数量

library(VennDiagram)
venn.diagram(list(sample_C = C,sample_D = D),fill = c("yellow","cyan"), cex = 1.5,
             filename = "venn2.png")
#    sampleC也可以是其他，fill是填充颜色，col是边框颜色，cex字体大小；
#    只能通过文件查看图像，无法输出到屏幕

dir()
#   查看当前目录下是否生成"venn2.png"

venn.diagram(list(A = A, C = C, D = D), fill = c("yellow","red","cyan"),
             cex = 1.5,filename="venn3.png")
venn.diagram(list(A = A, B = B, C = C, D = D), 
             fill = c("yellow","red","cyan","forestgreen"),
             cex = 1.5,filename="venn4.png")
venn.diagram(list(A = A, B = B, C = C, D = D , E = E ),
             fill = c("yellow","red","cyan","forestgreen","lightblue"),
             cex = 1.5,filename="venn5.png")
#    分别是2、3、4、5个元素/样品的韦恩图


#####################################
##        第15节 曼哈顿图          ##
#####################################

install.packages("qqman")
library(qqman)
library(RColorBrewer)
str(gwasResults)
#    查看数据类型

head(gwasResults)
#     查看数据，包括SNP名称，所处染色体编号，染色体坐标但因为是模拟数据所以瞎拍
#      P值

manhattan(gwasResults)
#    此时的图有两个阈值5,7

manhattan(gwasResults,main="Manhattan plot",ylim=c(0,6),cex=0.6)
#   设置标题，y轴大小，点的大小

manhattan(gwasResults, main = "Manhattan Plot", ylim = c(0, 6), cex = 0.6,
          cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline =F,
          genomewideline = F,chrlabs = c(1:20, "P", "Q"))

#    这些参数分别是x轴字体大小，圆点颜色，阈值线，染色体名称

manhattan(gwasResults,main="Manhttan plot",ylim=c(0,6),cex=0.6,cex.axis=0.9,
          col = c("blue4","orange3"),suggestiveline = 4.5,genomewideline = 5)
#    下面修改颜色

library(RColorBrewer)
color <- brewer.pal(name = "Set1",n = 4)
manhattan(gwasResults,main="Manhttan plot",ylim=c(0,6),cex=0.6,cex.axis=0.9,
          col = color,suggestiveline = 4.5,genomewideline = F)
#    可以看到超过阈值的点都在三号染色体，接下来对三号染色体进行单独分析

snp_3 <- subset(gwasResults,CHR==3)
head(snp_3)
#    取出3号染色体上的SNP

manhattan(snp_3)
#    可以看到300-400处分布明显

manhattan(subset(gwasResults,CHR==3))

#高亮显示部分SNP结果
snpsOfInterest
manhattan(gwasResults, highlight = snpsOfInterest)

#以下为注释SNP结果

manhattan(gwasResults, annotatePval = 0.001)
#   每条染色体上P值小于0.001的点中最明显的一个SNP

manhattan(gwasResults, annotatePval = 0.001, annotateTop = FALSE)
#   在上一个基础上所有满足条件的SNP

#更多内容可以查看manhattan与qqman的帮助文档
help("manhattan")
vignette("qqman")

#####################################
##        第16节 火山图            ##
#####################################

m <- read.csv("Deseq2.csv",header = T,row.names = 1)
head(m)
m <- na.omit(m)
#    去除NA缺失值

plot(m$log2FoldChange,m$padj)
plot(m$log2FoldChange,-1*log10(m$padj))
#    取padj-log10的值

plot(m$log2FoldChange,-1*log10(m$padj),xlim = c(-10,10),ylim=c(0,100))
#    设置x,y轴大小

#   接下来是拆分数据。火山图由三部分组成：上调表达基因、下调表达基因、不变基因
#   分别设置不同的颜色


m <- transform(m,padj=-1*log10(m$padj))
head(m)

#    小于-1为下调表达，大于1上调表达，二者之间表达量不变

down <- m[m$log2FoldChange<=-1,] 
up <- m[m$log2FoldChange>=1,]
no <- m[m$log2FoldChange>-1 & m$log2FoldChange <1,] 

plot(no$log2FoldChange,no$padj,xlim = c(-10,10),ylim=c(0,100),col="blue",
     pch=16,cex=0.8,main = "Gene Expression",
     xlab = "log2FoldChange",ylab="-log10(Qvalue)")

#    绘图，三部分依次加上去。每部分的log2FoldChange为x，padj为y。
#    pch=16表示实心点，cex点的大小

points(up$log2FoldChange,up$padj,col="red",pch=16,cex=0.8)
points(down$log2FoldChange,down$padj,col="green",pch=16,cex=0.8)


#####################################
##        第17节 GC-depth图        ##
#####################################

opar <- par(no.readonly = T)
par(mfrow=c(2,3))
#   绘制的多张图按2行3列排布

plot(pressure)
#   画六张，观察分布。它是平均分布，但这里的图不是，所以用layout设置

m <- read.table("GC-depth.txt");
head(m)
#   GC含量和测序深度

plot(m$V1,m$V2)
hist(m$V1)
hist(m$V2)
#   绘制对应的三部分，GC分布点图，GC直方图，深度直方图

nf <- layout(matrix(c(0,2,0,0,1,3),2,3,byrow=T),widths =  c(0.5,3,1),
             heights =  c(1,3,0.5),TRUE);
#    matrix表示在一个2行3列的格子内图片的布局，0表示无，1/2/3表示第1/2/3张图。
#    再设置每张图片的长宽比。最后为是否保留边距。

layout.show(3)
#    展示图片布局

par("mar")
#    查看边距

par(mar=c(5,5,0.5,0.5))
#    设置边距,分别代表下、左、上、右

layout.show(3)

x <- m[,1]
#    设置x轴

y <- m[,2]
#    设置y轴

plot(x,y,xlab='GC Content(%)',ylab='Depth',pch=16,col="#FF000077",
     xlim=c(0,100),ylim=c(0,max(y)))
#    先画散点图，即图1

hist(x,breaks = 100)
#    会显示x轴，影响美观

xhist <- hist(x,breaks=100,plot=FALSE);
yhist <- hist(y,breaks=floor(max(y)-0),plot=FALSE);
#    也可以先将图保存到变量里

xhist
xhist$counts
#    获得不同GC含量分布频数

#   此时需要重新布局

hist(x,breaks = 100)


plot(x,y,xlab='GC Content(%)',ylab='Depth',pch=16,col="#FF000077",
     xlim=c(0,100),ylim=c(0,max(y)))
#    重新绘制散点图

par(mar=c(0,5,1,1))
#    对2布局

barplot(xhist$counts,space=0,xlim=c(0,100) )

par(mar=c(5,0,1,1))
#   对3布局

barplot(yhist$counts,space=0,horiz=TRUE,ylim=c(0,max(y) ) );


#####################################
##        第18节 COG注释图         ##
#####################################

m <- read.table("cog.class.annot.txt",head=T,sep="\t")
head(m)
#    由3列组成，第一列为大类，第二列是具体功能，第三列是数量

layout(matrix(c(1,2),nrow = 1),widths=c(20,13))
#    布局图像，1行两列，宽为20：13

layout.show(2)
par( mar=c(3,4,4,1)+0.1 )
#    展示布局

class <- c("J","A","K","L","B","D","Y","V","T","M","N","Z","W"
           ,"U","O","C","G","E","F","H","I","P","Q","R","S")
class
#    根据数据框比对结果将每个大类排序，系统默认按ABC排序

factor( as.character(m$Code),levels=class )
t <- factor( as.character(m$Code),levels=class )
t
#    根据大类转换为因子，class排序结果作为水平

m <- m[order(t),]
head(m)
#    重新对数据进行排序

barplot(m$Gene.Number,space=F,col=rainbow(25),
        ylab="Number of genes",names.arg = m$Code)
#    space条形图之间有无空隙


#    接下来根据比对结果添加注释

l <- c(0,5,15,23,25)
#    1-5，6-15，16-23，24-25共4类

id<- c("INFORMATION STORAGE\nAND PROCESSING",
       "CELLULAR PROCESSES\nAND SIGNALING",
       "METABOLISM","POORLY\nCHARACTERIZED")

#    根据结果添加每一类型的名称

abline( v = l[c(-1,-5)])
#    添加分组线，v表示垂直线,后面表示垂线添加的位置，除了x=0/25外，其它3处添加

for( i in 2:length(l) ){
  text( (l[i-1]+l[i])/2,max(m[,3])*1.1,id[i-1],cex=0.8,xpd=T)
}

#    循环为2-5即i=2-5，每次横坐标在所处区域正中，纵坐标为第三行数据最大值的1.1倍
#     每次写第i-1个tittle，字号为0.8，是否显示标题。

par(mar=c(2,0,2,1)+0.1 )
#    设置图2边距，下左上右各加0.1

plot(0,0,type="n",xlim=c(0,1),ylim=c(0,26),bty="n",axes=F,xlab="",ylab="")
#    绘制一个空图，x/y大小，aex=F表示不显示x/y轴。

for( i in 1:length(class) ){
  text(0,26-i+0.5,paste(m$Code[i],m$Functional.Categories[i]),pos=4,cex=1,pty=T)
}
#    text每次在横坐标为0，纵坐标为26-i.5处写入第i个大类和功能
#    pos上中等于1，右中为2，下中为3，左中为4

title(main = "COG function classification");
#    图标题



#####################################
##        第19节 GO条目图          ##
#####################################

library(ggplot2)

go <- read.csv("Rdata/go.csv",header = T)
head(go)
View(go)
#    这里head不方便观察，因此换做view，选1，2，5列绘图

go_sort <- go[order(go$Ontology,-go$Percentage),]
#    order排序，-表降序排列

m <- go_sort[go_sort$Ontology=="Molecular function",][1:10,]
c <- go_sort[go_sort$Ontology=="Cellular component",][1:10,]
b <- go_sort[go_sort$Ontology=="Biological process",][1:10,]
#    GO注释有Molecular Function、Cellular component、biology process 3个term
#    每个term有n个数据，只取前10-20个

slimgo <- rbind(b,c,m)
#    rbind表按行合并

slimgo$Term=factor(slimgo$Term,levels=slimgo$Term)
#    首先需要将Trem转换为因子

colnames(slimgo)
pp <- ggplot(data = slimgo, mapping = aes(x=Term,y=Percentage,fill=Ontology))
#    ggplot绘图，数据集，Term为x，percentage为y，Ontology作为填充色的分类

pp
pp+geom_bar(stat="identity")
#    绘制分组直方图

pp+geom_bar(stat="identity")+coord_flip()
#    翻转坐标轴

pp+geom_bar(stat="identity")+coord_flip()+
  scale_x_discrete(limits=rev(levels(slimgo$Term)))
#    scale_x_discrete 表x不是连续变量，rev表反转排序

pp+geom_bar(stat="identity")+coord_flip()+
  scale_x_discrete(limits=rev(levels(slimgo$Term)))+guides(fill=FALSE)
#    不显示图例，会出现警告，但不是错误

pp+geom_bar(stat="identity")+coord_flip()+
  scale_x_discrete(limits=rev(levels(slimgo$Term)))+
  guides(fill=FALSE)+theme_bw()
#     将背景色有默认的灰色改为黑白色

#####################################
##        第20节 kegg图            ##
#####################################

library(ggplot2)

pathway <-  read.csv("Rdata/kegg.csv",header=T)
View(pathway)
#    Pathway通路名称、AvsB差异表达基因数量、All_Unigene参与该通路的所有基因
#    Pvalue P值、Qvalue Q值、richFactor=All_Unigene/AvsB、Pathway.ID 通路ID
#    Genes  通路中的基因ID、KOs 通路中的基因qq号？
#    前6个最重要！！！

colnames(pathway)

pp <-  ggplot(data=pathway,mapping = aes(x=richFactor,y=Pathway))
#    pathway数据集，richFactor为x轴，Pathway为y轴

pp
pp + geom_point()
#    点图
pp + geom_point(aes(size=AvsB))
#    点的大小为AvsB差异表达基因数量

pp + geom_point(aes(size=AvsB,color=Qvalue))
#    颜色根据Q值变化而变化

pp + geom_point(aes(size=AvsB,color=Qvalue)) + 
  scale_colour_gradient(low="green",high="red")
#    改变颜色

pp + geom_point(aes(size=AvsB,color=Qvalue)) +
  scale_colour_gradient(low="green",high="red")+
  labs(title="Top20 of pathway enrichment",x="Rich factor",
       y="Pathway Name",color="-log10(Qvalue)",size="Gene Numbers")

#    添加标签

pp + geom_point(aes(size=AvsB,color=Qvalue)) + 
  scale_colour_gradient(low="green",high="red")+
  labs(title="Top20 of pathway enrichment",x="Rich factor",y="Pathway Name",
       color="-log10(Qvalue)",size="Gene Numbers")+theme_bw()
#    改变背景色
