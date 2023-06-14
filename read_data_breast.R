if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("curatedBreastData")
library(curatedBreastData)


# read the data, study id 2034

data(curatedBreastDataExprSetList)
D=curatedBreastDataExprSetList$study_2034_GPL96_all

Data=data.frame(Y=D$RFS_months_or_MIN_months_of_RFS,status=1-D$RFS,age=D$age,ER=D$ER_preTrt,sizeTum=D$tumor_stage_preTrt,menopause=D$menopausal_status)
Data$sizeTum[which(Data$sizeTum=='T1')]=1
Data$sizeTum[which(Data$sizeTum=='T2')]=2
Data$sizeTum[which(Data$sizeTum=='T3')]=3
Data$sizeTum[which(Data$sizeTum=='T4')]=4

# remove those with Tumour size > 2 (in total, 8 removed)
Data_breast=Data[-which(Data$sizeTum>2),]

# standardize age
#age_backup <- Data_breast$age
Data_breast$age <- Data_breast$age - mean(Data_breast$age)


n=dim(Data_breast)[1]    # sample size
Y=Data_breast[,1]        # follow-up time
status=Data_breast[,2]   # censoring status:  censored (status=0), noncensored (status=1)

X=Data_breast[,c(3:6)]
Z=Data_breast[,c(3:6)]   

# ER+
Data1=data.frame(Y=Y[X[,2]==1], status=status[X[,2]==1], age=Data_breast$age[X[,2]==1], sizeTum=Data_breast$sizeTum[X[,2]==1])
# ER-
Data2=data.frame(Y=Y[X[,2]==0], status=status[X[,2]==0], age=Data_breast$age[X[,2]==0], sizeTum=Data_breast$sizeTum[X[,2]==0])#, post_men=X$post_men[X[,2]==0])

