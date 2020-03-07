########## Load library #####
library(dplyr)
library(GRS.test)

########## Define start and end date for calculation #########
start.date = 196307
end.date = 201412
column.count = 618

## Define functions for data processing
process.input.data <- function(input.df) {
  colnames(input.df)[1] <- "YYYYMM"
  n.col = ncol(input.df)
  for (i in 1:n.col) {
    input.df[,i] <- as.numeric(as.character(input.df[,i]))}
  return(input.df)
}

filter.by.date <- function(input.df) {
  output.df <- NULL
  output.df <- filter(input.df,as.numeric(YYYYMM)>=start.date & as.numeric(YYYYMM)<=end.date)
  return(output.df)
}

########################################################### Table 1 ###################################################################

# Read and process factors: RM-RF SMB HML RMW CMA
five.factors <- process.input.data(read.csv("F-F_Research_Data_5_Factors.CSV",header=TRUE))
five.factors.filtered <- filter.by.date(five.factors)

# Read and process factor MOM
factor.mom <- process_input_data(read.csv("F-F_Momentum_Factor.CSV",header=TRUE))
factor.mom.filtered <- filter.by.date(factor.mom)

# Formulate Table 1
five.factors.filtered$MOM <- factor.mom.filtered$Mom
TABLE.ONE <- rbind (round(apply(five.factors.filtered[c(2:6, 8)], 2, mean, na.rm=TRUE), digits=2),
                    round(apply(five.factors.filtered[c(2:6,8)], 2, sd, na.rm=TRUE), digits=2),
                    round(apply(five.factors.filtered[c(2:6,8)], 2, function(x) t.test(x,mu=0)$statistic), digits=2)
                    )
rownames(TABLE.ONE) = c('Mean', 'SD', 't-statistic')
print(TABLE.ONE)
#write.csv(TABLE.ONE, "T1.CSV")

########################################################### Table 2 ###################################################################

# Select Risk Free Return
RF <- five.factors.filtered$RF

# Define colume name
TABLE.TWO.COL.NAMES <- c("GRS",
                         "p(GRS)",
                         "A|ai|",
                         "A|ai|/A|ri|",
                         "A(ai^2)/A(ri^2)",
                         "A(s^2(ai)/A(ri^2)",
                         "A(R^2)")

# Function to generate one row for a model
generate.t2.row <- function(input.ret, factor.string) {
  return.mat <- input.ret
  no.of.factor <- length(factor.string)
  factor.mat = data.matrix(five.factors.filtered[, factor.string])
  result <- GRS.test(return.mat, factor.mat)
  
  # calculate alpha
  reg.func <- function (y, m) {
    mod <- lm(y ~ m)
    mod$coefficients["(Intercept)"]
  }
  alphas <- as.vector(apply(factor.mat, 2, reg.func, return.mat))
  
  # calculate A(ai^2)/A(ri^2)
  Ri <- colMeans(return.mat) - mean(five.factors.filtered$Mkt.RF)
  reg.func.error <- function (y, m) {
    mod <- lm(y ~ m)
    summary(mod)$coefficients[1,2]
  }
  alphas.error <- as.vector(apply(return.mat, 2, reg.func.error, factor.mat))
  
  # calculate adjusted R squared
  reg.func.adj.R.squared <- function (y, m) {
    mod <- lm(y ~ m)
    summary(mod)$adj.r.squared
  }
  alphas.adj.R.squared <- as.vector(apply(return.mat, 2, reg.func.adj.R.squared, factor.mat))
  
  output.row <- c(result$GRS.stat,
                  result$GRS.pval,
                  mean(abs(alphas)),
                  mean(abs(alphas))/mean(abs(Ri)),
                  mean(alphas^2)/mean(Ri^2),
                  mean(alphas.error^2)/mean(alphas^2),
                  mean(alphas.adj.R.squared)
  )
}

###################################
##### 25 Size-beta portfolios #####
size.beta <- read.csv("25_Portfolios_ME_BETA.CSV",header=TRUE)
size.beta <- size.beta[1:column.count,] 
size.beta <- process.input.data(size.beta)
size.beta.filtered <- filter.by.date(size.beta)

# Form GRS test factor and return matrix
return.mat.beta <- data.matrix(size.beta.filtered[,2:26])
return.mat.beta <- sweep(return.mat.beta, 2, RF, FUN='-')

size.beta.Mkt <- generate.t2.row(return.mat.beta, c("Mkt.RF")) 
size.beta.Mkt.SMB.HML <- generate.t2.row(return.mat.beta, c("Mkt.RF", "SMB", "HML")) 
size.beta.Mkt.SMB.HML.RMW <- generate.t2.row(return.mat.beta, c("Mkt.RF", "SMB", "HML", "RMW")) 
size.beta.Mkt.SMB.HML.CMA <- generate.t2.row(return.mat.beta, c("Mkt.RF", "SMB", "HML", "CMA")) 
size.beta.Mkt.SMB.RMW.CMA <- generate.t2.row(return.mat.beta, c("Mkt.RF", "SMB", "RMW", "CMA")) 
size.beta.Mkt.SMB.HML.RMW.CMA <- generate.t2.row(return.mat.beta, c("Mkt.RF", "SMB", "HML", "RMW", "CMA")) 

#################################
##### 35 Size-NI portfolios #####
size.NI <- read.csv("35_Portfolios_ME_NI.CSV",header=TRUE)
size.NI <- size.NI[1:column.count,] 
size.NI <- process.input.data(size.NI)
size.NI.filtered <- filter.by.date(size.NI)

# Form GRS test factor and return matrix
return.mat.NI <- data.matrix(size.NI.filtered[,2:36])
return.mat.NI <- sweep(return.mat.NI, 2, RF, FUN='-')

size.NI.Mkt.SMB.HML <- generate.t2.row(return.mat.NI, c("Mkt.RF", "SMB", "HML"))
size.NI.Mkt.SMB.HML.RMW <- generate.t2.row(return.mat.NI, c("Mkt.RF", "SMB", "HML", "RMW"))
size.NI.Mkt.SMB.HML.CMA <- generate.t2.row(return.mat.NI, c("Mkt.RF", "SMB", "HML", "CMA"))
size.NI.Mkt.SMB.RMW.CMA <- generate.t2.row(return.mat.NI, c("Mkt.RF", "SMB", "RMW", "CMA"))
size.NI.Mkt.SMB.HML.RMW.CMA <- generate.t2.row(return.mat.NI, c("Mkt.RF", "SMB", "HML", "RMW", "CMA"))

##################################
##### 25 Size-Var portfolios #####
size.var <- read.csv("25_Portfolios_ME_VAR.CSV",header=TRUE)
size.var <- size.var[1:column.count,] 
size.var <- process.input.data(size.var)
size.var.filtered <- filter.by.date(size.var)

# Form GRS test factor and return matrix
return.mat.var <- data.matrix(size.var.filtered[,2:26])
return.mat.var <- sweep(return.mat.var, 2, RF, FUN='-')

size.var.Mkt.SMB.HML <- generate.t2.row(return.mat.var, c("Mkt.RF", "SMB", "HML"))
size.var.Mkt.SMB.HML.RMW <- generate.t2.row(return.mat.var, c("Mkt.RF" , "SMB", "HML", "RMW"))
size.var.Mkt.SMB.HML.CMA <- generate.t2.row(return.mat.var, c("Mkt.RF" , "SMB", "HML", "CMA"))
size.var.Mkt.SMB.RMW.CMA <- generate.t2.row(return.mat.var, c("Mkt.RF" , "SMB", "RMW", "CMA"))
size.var.Mkt.SMB.HML.RMW.CMA <- generate.t2.row(return.mat.var, c("Mkt.RF", "SMB", "HML", "RMW", "CMA"))

###########################################
##### 25 Size-Residual Var portfolios #####
size.Rvar <- read.csv("25_Portfolios_ME_RESVAR.CSV", header=TRUE)
size.Rvar <- size.Rvar[1:column.count,] 
size.Rvar <- process.input.data(size.Rvar)
size.Rvar.filtered <- filter.by.date(size.Rvar)

# Form GRS test factor and return matrix
return.mat.Rvar <- data.matrix(size.Rvar.filtered[,2:26])
return.mat.Rvar <- sweep(return.mat.Rvar, 2, RF, FUN='-')

size.Rvar.Mkt.SMB.HML <- generate.t2.row(return.mat.Rvar, c("Mkt.RF", "SMB", "HML"))
size.Rvar.Mkt.SMB.HML.RMW <- generate.t2.row(return.mat.Rvar, c("Mkt.RF", "SMB", "HML", "RMW"))
size.Rvar.Mkt.SMB.HML.CMA <- generate.t2.row(return.mat.Rvar, c("Mkt.RF", "SMB", "HML", "CMA"))
size.Rvar.Mkt.SMB.RMW.CMA <- generate.t2.row(return.mat.Rvar, c("Mkt.RF", "SMB", "RMW", "CMA"))
size.Rvar.Mkt.SMB.HML.RMW.CMA <- generate.t2.row(return.mat.Rvar, c("Mkt.RF", "SMB", "HML", "RMW", "CMA"))

#################################
##### 25 Size-AC portfolios #####
size.AC <- read.csv("25_Portfolios_ME_AC.CSV",header=TRUE)
size.AC <- size.AC[1:column.count,] 
size.AC <- process.input.data(size.AC)
size.AC.filtered <- filter.by.date(size.AC)

# Form GRS test factor and return matrix
return.mat.AC <- data.matrix(size.AC.filtered[,2:26])
return.mat.AC <- sweep(return.mat.AC, 2, RF, FUN='-')

size.AC.Mkt.SMB.HML <- generate.t2.row(return.mat.AC, c("Mkt.RF", "SMB", "HML"))
size.AC.Mkt.SMB.HML.RMW <- generate.t2.row(return.mat.AC, c("Mkt.RF", "SMB", "HML", "RMW"))
size.AC.Mkt.SMB.HML.CMA <- generate.t2.row(return.mat.AC, c("Mkt.RF", "SMB", "HML", "CMA"))
size.AC.Mkt.SMB.RMW.CMA <- generate.t2.row(return.mat.AC, c("Mkt.RF", "SMB", "RMW", "CMA"))
size.AC.Mkt.SMB.HML.RMW.CMA <- generate.t2.row(return.mat.AC, c("Mkt.RF", "SMB", "HML", "RMW", "CMA"))

#########################################
##### 25 Size-Prior 2-12 portfolios #####
size.prior <- read.csv("25_Portfolios_ME_MOM.CSV",header=TRUE)
size.prior <- size.prior[1:1106,] 
size.prior <- process.input.data(size.prior)
size.prior.filtered <- filter.by.date(size.prior)

# Form GRS test factor and return matrix
return.mat.prior <- data.matrix(size.prior.filtered[,2:26])
return.mat.prior <- sweep(return.mat.prior, 2, RF, FUN='-')

size.prior.Mkt.SMB.HML <- generate.t2.row(return.mat.prior, c("Mkt.RF", "SMB", "HML"))
size.prior.Mkt.SMB.HML.RMW <- generate.t2.row(return.mat.prior, c("Mkt.RF", "SMB", "HML", "RMW"))
size.prior.Mkt.SMB.HML.CMA <- generate.t2.row(return.mat.prior, c("Mkt.RF", "SMB", "HML", "CMA"))
size.prior.Mkt.SMB.RMW.CMA <- generate.t2.row(return.mat.prior, c("Mkt.RF", "SMB", "RMW" , "CMA"))
size.prior.Mkt.SMB.HML.RMW.CMA <- generate.t2.row(return.mat.prior, c("Mkt.RF", "SMB", "HML", "RMW", "CMA"))

size.prior.Mkt.SMB.HML.MOM <- generate.t2.row(return.mat.prior, c("Mkt.RF", "SMB", "HML", "MOM"))
size.prior.Mkt.SMB.HML.RMW.MOM <- generate.t2.row(return.mat.prior, c("Mkt.RF", "SMB", "HML", "RMW", "MOM"))
size.prior.Mkt.SMB.HML.CMA.MOM <- generate.t2.row(return.mat.prior, c("Mkt.RF", "SMB", "HML", "CMA", "MOM"))
size.prior.Mkt.SMB.RMW.CMA.MOM <- generate.t2.row(return.mat.prior, c("Mkt.RF", "SMB", "RMW", "CMA", "MOM"))
size.prior.Mkt.SMB.HML.RMW.CMA.MOM <- generate.t2.row(return.mat.prior, c("Mkt.RF", "SMB", "HML", "RMW", "CMA", "MOM"))

##############################
TABLE.TWO <-data.frame(rbind(
  size.beta.Mkt,
  size.beta.Mkt.SMB.HML,
  size.beta.Mkt.SMB.HML.RMW,
  size.beta.Mkt.SMB.HML.CMA,
  size.beta.Mkt.SMB.RMW.CMA,
  size.beta.Mkt.SMB.HML.RMW.CMA,
  
  size.NI.Mkt.SMB.HML,
  size.NI.Mkt.SMB.HML.RMW, 
  size.NI.Mkt.SMB.HML.CMA,
  size.NI.Mkt.SMB.RMW.CMA,
  size.NI.Mkt.SMB.HML.RMW.CMA,
  
  size.var.Mkt.SMB.HML,
  size.var.Mkt.SMB.HML.RMW,
  size.var.Mkt.SMB.HML.CMA,
  size.var.Mkt.SMB.RMW.CMA,
  size.var.Mkt.SMB.HML.RMW.CMA,
  
  size.Rvar.Mkt.SMB.HML,
  size.Rvar.Mkt.SMB.HML.RMW,
  size.Rvar.Mkt.SMB.HML.CMA,
  size.Rvar.Mkt.SMB.RMW.CMA,
  size.Rvar.Mkt.SMB.HML.RMW.CMA,
  
  size.AC.Mkt.SMB.HML,
  size.AC.Mkt.SMB.HML.RMW,
  size.AC.Mkt.SMB.HML.CMA,
  size.AC.Mkt.SMB.RMW.CMA,
  size.AC.Mkt.SMB.HML.RMW.CMA,
  
  size.prior.Mkt.SMB.HML,
  size.prior.Mkt.SMB.HML.RMW ,
  size.prior.Mkt.SMB.HML.CMA ,
  size.prior.Mkt.SMB.RMW.CMA ,
  size.prior.Mkt.SMB.HML.RMW.CMA ,
  
  size.prior.Mkt.SMB.HML.MOM ,
  size.prior.Mkt.SMB.HML.RMW.MOM ,
  size.prior.Mkt.SMB.HML.CMA.MOM ,
  size.prior.Mkt.SMB.RMW.CMA.MOM ,
  size.prior.Mkt.SMB.HML.RMW.CMA.MOM))

colnames(TABLE.TWO) = TABLE.TWO.COL.NAMES
print(TABLE.TWO)

#write.csv(TABLE.TWO, "T2.CSV")

########################################################### Table 3 ###################################################################

###########################################
##### 25 Size-B/M portfolios #####

TABLE.THREE.COL.NAMES <- c("Mean: Low", "2", "3", "4", "High", "SD: Low", "2", "3", "4", "High")

size.bm <- read.csv("25_Portfolios_ME_BM.CSV",header=TRUE)
size.bm <- size.bm[1:column.count,] 
size.bm <- process.input.data(size.bm)
size.bm.filtered <- filter.by.date(size.bm)
size.bm.excess = sweep(size.bm.filtered[,2:26], 2, RF, FUN="-")

size.bm.mean = t(as.matrix(round(apply(size.bm.excess, 2, mean), 2)))
size.bm.sd = t(as.matrix(round(apply(size.bm.excess, 2, sd), 2)))

Small = t(as.matrix(c(size.bm.mean[1, 1:5], size.bm.sd[1, 1:5])))
Second = t(as.matrix(c(size.bm.mean[1, 6:10], size.bm.sd[1, 6:10])))
Third = t(as.matrix(c(size.bm.mean[1, 11:15], size.bm.sd[1, 11:15])))
Fourth = t(as.matrix(c(size.bm.mean[1, 16:20], size.bm.sd[1, 16:20])))
Fifth = t(as.matrix(c(size.bm.mean[1, 21:25], size.bm.sd[1, 21:25])))
colnames(Small) = TABLE.THREE.COL.NAMES
colnames(Second) = TABLE.THREE.COL.NAMES
colnames(Third) = TABLE.THREE.COL.NAMES
colnames(Fourth) = TABLE.THREE.COL.NAMES
colnames(Fifth) = TABLE.THREE.COL.NAMES


TABLE.THREE <- data.frame(rbind(
  Small, Second, Third, Fourth, Fifth
))

colnames(TABLE.THREE) = TABLE.THREE.COL.NAMES
rownames(TABLE.THREE) = c("Small", "2", "3", "4", "Big")
print(TABLE.THREE)

# write.csv(TABLE.THREE, "T3.CSV")

########################################################### Table 4 ###################################################################
regression.line <- function (y, m) {
  mod <- lm(y ~ m)
  summary(mod)$coefficients
}

size.beta.excess <- as.matrix(sweep(data.matrix(size.beta.filtered[,2:26]), 2, RF, FUN='-'))
MKT.RF = as.matrix(five.factors.filtered[, 2])
SMB = as.matrix(five.factors.filtered[, 3])
HML = as.matrix(five.factors.filtered[, 4])
RMW = as.matrix(five.factors.filtered[, 5])
CMA = as.matrix(five.factors.filtered[, 6])

T4.Panel.A <- matrix(-99, 10, 10)
for(i in 1 : 25) {
  summary = regression.multiple(size.beta.excess[, i], MKT.RF)
  T4.Panel.A[((i-1) %/% 5 + 1), ((i- 1) %% 5 + 1)] = round(summary[1, 1], 2)
  T4.Panel.A[((i-1) %/% 5 + 1), ((i- 1) %% 5 + 1) + 5] = round(summary[1, 3], 2)
  T4.Panel.A[((i-1) %/% 5 + 1) + 5, ((i- 1) %% 5 + 1)] = round(summary[2, 1], 2)
  T4.Panel.A[((i-1) %/% 5 + 1) + 5, ((i- 1) %% 5 + 1) + 5] = round(summary[2, 3], 2)
}
colnames(T4.Panel.A) <- c("low", "2", "3", "4", "high", "low", "2", "3", "4", "high")
rownames(T4.Panel.A) <- c("small", "2", "3", "4", "big", "small", "2", "3", "4", "big")
print(T4.Panel.A)

T4.Panel.B <- matrix(-99, 30, 10)
for(i in 1 : 25) {
  y = as.matrix(size.beta.excess[, i])
  y.x = y.x = cbind(y, five.factors.filtered[, 2:6])
  summary = summary(lm(formula = y ~ MKT.RF + SMB + HML + RMW + CMA, data = y.x))$coefficients
  T4.Panel.B[((i-1) %/% 5 + 1), ((i- 1) %% 5 + 1)] = round(summary[1, 1], 2)
  T4.Panel.B[((i-1) %/% 5 + 1) + 5, ((i- 1) %% 5 + 1)] = round(summary[2, 1], 2)
  T4.Panel.B[((i-1) %/% 5 + 1) + 10, ((i- 1) %% 5 + 1)] = round(summary[3, 1], 2)
  T4.Panel.B[((i-1) %/% 5 + 1) + 15, ((i- 1) %% 5 + 1)] = round(summary[4, 1], 2)
  T4.Panel.B[((i-1) %/% 5 + 1) + 20, ((i- 1) %% 5 + 1)] = round(summary[5, 1], 2)
  T4.Panel.B[((i-1) %/% 5 + 1) + 25, ((i- 1) %% 5 + 1)] = round(summary[6, 1], 2)
  
  T4.Panel.B[((i-1) %/% 5 + 1), ((i- 1) %% 5 + 1) + 5] = round(summary[1, 3], 2)
  T4.Panel.B[((i-1) %/% 5 + 1) + 5, ((i- 1) %% 5 + 1) + 5] = round(summary[2, 3], 2)
  T4.Panel.B[((i-1) %/% 5 + 1) + 10, ((i- 1) %% 5 + 1) + 5] = round(summary[3, 3], 2)
  T4.Panel.B[((i-1) %/% 5 + 1) + 15, ((i- 1) %% 5 + 1) + 5] = round(summary[4, 3], 2)
  T4.Panel.B[((i-1) %/% 5 + 1) + 20, ((i- 1) %% 5 + 1) + 5] = round(summary[5, 3], 2)
  T4.Panel.B[((i-1) %/% 5 + 1) + 25, ((i- 1) %% 5 + 1) + 5] = round(summary[6, 3], 2)
  
}

colnames(T4.Panel.B) <- c("low", "2", "3", "4", "high", "low", "2", "3", "4", "high")
rownames(T4.Panel.B) <- c("small", "2", "3", "4", "big", "small", "2", "3", "4", "big", 
                          "small", "2", "3", "4", "big", "small", "2", "3", "4", "big", 
                          "small", "2", "3", "4", "big", "small", "2", "3", "4", "big")

print(T4.Panel.B)
TABLE.FOUR = rbind(
  c("PanelA: CAPM: Rit – RFt = ai + bi(RMt – RFt) + eit ", "", "", "", "", "", "", "", "", ""),
  c("", "", "a", "", "", "", "", "t(a)", "", ""), 
  T4.Panel.A[1:5,], 
  c("", "", "b", "", "", "", "", "t(b)", "", ""),
  T4.Panel.A[6:10,],
  c("Panel B: Five-factor: Rit – RFt = ai + bi(RMt – RFt) + siSMBt + hiHMLOt + riRMWt + ciCMAt + eit ", "", "", "", "", "", "", "", "", ""),
  c("", "", "a", "", "", "", "", "t(a)", "", ""), 
  T4.Panel.B[1:5,], 
  c("", "", "b", "", "", "", "", "t(b)", "", ""),
  T4.Panel.B[6:10,], 
  c("", "", "s", "", "", "", "", "t(s)", "", ""), 
  T4.Panel.B[11:15,], 
  c("", "", "h", "", "", "", "", "t(h)", "", ""), 
  T4.Panel.B[16:20,], 
  c("", "", "r", "", "", "", "", "t(r)", "", ""), 
  T4.Panel.B[21:25,], 
  c("", "", "c", "", "", "", "", "t(c)", "", ""),
  T4.Panel.B[26:30,]
  )

#write.csv(TABLE.FOUR, "T4.CSV")

########################################################### Table 5 ###################################################################

###########################################
##### 35 Size-NI portfolios #####

TABLE.FIVE.COL.NAMES <- c("Low", "2", "3", "4", "High", "Low", "2", "3", "4", "High")

size.ni <- read.csv("35_Portfolios_ME_NI.CSV",header=TRUE)
size.ni <- size.ni[1:column.count,]
size.ni <- process.input.data(size.ni)
size.ni.filtered <- filter.by.date(size.ni)
size.ni.excess = sweep(size.ni.filtered[,2:26], 2, RF, FUN="-")

size.ni.mean = t(as.matrix(round(apply(size.ni.excess, 2, mean), 2)))
size.ni.sd = t(as.matrix(round(apply(size.ni.excess, 2, sd), 2)))

Small = t(as.matrix(c(size.ni.mean[1, 1:5], size.ni.sd[1, 1:5])))
Second = t(as.matrix(c(size.ni.mean[1, 6:10], size.ni.sd[1, 6:10])))
Third = t(as.matrix(c(size.ni.mean[1, 11:15], size.ni.sd[1, 11:15])))
Fourth = t(as.matrix(c(size.ni.mean[1, 16:20], size.ni.sd[1, 16:20])))
Fifth = t(as.matrix(c(size.ni.mean[1, 21:25], size.ni.sd[1, 21:25])))
colnames(Small) = TABLE.FIVE.COL.NAMES
colnames(Second) = TABLE.FIVE.COL.NAMES
colnames(Third) = TABLE.FIVE.COL.NAMES
colnames(Fourth) = TABLE.FIVE.COL.NAMES
colnames(Fifth) = TABLE.FIVE.COL.NAMES


TABLE.FIVE <- data.frame(rbind(
  Small, Second, Third, Fourth, Fifth
))

colnames(TABLE.FIVE) = TABLE.FIVE.COL.NAMES
rownames(TABLE.FIVE) = c("Small", "2", "3", "4", "Big")
print(TABLE.FIVE)
#write.csv(TABLE.FIVE, "T5.CSV")

########################################################### Table 6 ###################################################################
size.ni.excess <- as.matrix(sweep(data.matrix(size.ni.filtered[,2:36]), 2, RF, FUN='-'))
MKT.RF = as.matrix(five.factors.filtered[, 2])
SMB = as.matrix(five.factors.filtered[, 3])
HML = as.matrix(five.factors.filtered[, 4])
RMW = as.matrix(five.factors.filtered[, 5])
CMA = as.matrix(five.factors.filtered[, 6])

T6.Panel.A <- matrix(-99, 5, 14)
for(i in 1 : 35) {
  y = as.matrix(size.ni.excess[, i])
  y.x = cbind(y, five.factors.filtered[, 2:6])
  summary = summary(lm(formula = y ~ MKT.RF + SMB + HML, data = y.x))$coefficients
  T6.Panel.A[((i-1) %/% 7 + 1), ((i- 1) %% 7 + 1)] = round(summary[1, 1], 2)
  T6.Panel.A[((i-1) %/% 7 + 1), ((i- 1) %% 7 + 1) + 7] = round(summary[1, 3], 2)
}
colnames(T6.Panel.A) <- c("neg", "zero", "low", "2", "3", "4", "high", "neg", "zero", "low", "2", "3", "4", "high")
rownames(T6.Panel.A) <- c("small", "2", "3", "4", "big")
print(T6.Panel.A)

T6.Panel.B <- matrix(-99, 20, 14)
for(i in 1 : 35) {
  y = as.matrix(size.ni.excess[, i])
  y.x = cbind(y, five.factors.filtered[, 2:6])
  summary = summary(lm(formula = y ~ MKT.RF + SMB + HML + RMW + CMA, data = y.x))$coefficients
  T6.Panel.B[((i-1) %/% 7 + 1), ((i- 1) %% 7 + 1)] = round(summary[1, 1], 2)
  T6.Panel.B[((i-1) %/% 7 + 1) + 5, ((i- 1) %% 7 + 1)] = round(summary[4, 1], 2)
  T6.Panel.B[((i-1) %/% 7 + 1) + 10, ((i- 1) %% 7 + 1)] = round(summary[5, 1], 2)
  T6.Panel.B[((i-1) %/% 7 + 1) + 15, ((i- 1) %% 7 + 1)] = round(summary[6, 1], 2)
  
  T6.Panel.B[((i-1) %/% 7 + 1), ((i- 1) %% 7 + 1) + 7] = round(summary[1, 3], 2)
  T6.Panel.B[((i-1) %/% 7 + 1) + 5, ((i- 1) %% 7 + 1) + 7] = round(summary[4, 3], 2)
  T6.Panel.B[((i-1) %/% 7 + 1) + 10, ((i- 1) %% 7 + 1) + 7] = round(summary[5, 3], 2)
  T6.Panel.B[((i-1) %/% 7 + 1) + 15, ((i- 1) %% 7 + 1) + 7] = round(summary[6, 3], 2)
}

colnames(T6.Panel.B) <- c("neg", "zero", "low", "2", "3", "4", "high", "neg", "zero", "low", "2", "3", "4", "high")
rownames(T6.Panel.B) <- c("small", "2", "3", "4", "big", "small", "2", "3", "4", "big", 
                          "small", "2", "3", "4", "big", "small", "2", "3", "4", "big")

print(T6.Panel.B)

# rbind(T6.Panel.A, T6.Panel.B)

TABLE.SIX = rbind(
  c("PanelA:Three-factor:Rit–RFt=ai+bi(RMt–RFt)+siSMBt+hiHMLt+eit", "", "", "", "", "", "", "", "", "", "", "", "", ""),
  c("", "", "", "a", "", "", "", "", "", "", "t(a)", "", "", ""), 
  T6.Panel.A[1:5,], 
  c("Panel B: Five-factor: Rit – RFt = ai + bi(RMt – RFt) + siSMBt + hiHMLOt + riRMWt + ciCMAt + eit ", "", "", "", "", "", "", "", "", "", "", "", "", ""),
  c("", "", "", "a", "", "", "", "", "", "", "t(a)", "", "", ""),
  T6.Panel.B[1:5,],
  c("", "", "", "h", "", "", "", "", "", "", "t(h)", "", "", ""), 
  T6.Panel.B[6:10,], 
  c("", "", "", "r", "", "", "", "", "", "", "t(r)", "", "", ""),
  T6.Panel.B[11:15,],
  c("", "", "", "c", "", "", "", "", "", "", "t(c)", "", "", ""),
  T6.Panel.B[16:20,]
)
#write.csv(TABLE.SIX, "T6.CSV")

########################################################### Table 7 ###################################################################

###########################################
##### 25 Size-Rvar portfolios #####

TABLE.SEVEN.COL.NAMES <- c("Low", "2", "3", "4", "High", "Low", "2", "3", "4", "High")

size.rvar <- read.csv("25_Portfolios_ME_RESVAR.CSV",header=TRUE)
size.rvar <- size.rvar[1:column.count,] 
size.rvar <- process.input.data(size.rvar)
size.rvar.filtered <- filter.by.date(size.rvar)
size.rvar.excess = sweep(size.rvar.filtered[,2:26], 2, RF, FUN="-")

size.rvar.mean = t(as.matrix(round(apply(size.rvar.excess, 2, mean), 2)))
size.rvar.sd = t(as.matrix(round(apply(size.rvar.excess, 2, sd), 2)))

Small = t(as.matrix(c(size.rvar.mean[1, 1:5], size.rvar.sd[1, 1:5])))
Second = t(as.matrix(c(size.rvar.mean[1, 6:10], size.rvar.sd[1, 6:10])))
Third = t(as.matrix(c(size.rvar.mean[1, 11:15], size.rvar.sd[1, 11:15])))
Fourth = t(as.matrix(c(size.rvar.mean[1, 16:20], size.rvar.sd[1, 16:20])))
Fifth = t(as.matrix(c(size.rvar.mean[1, 21:25], size.rvar.sd[1, 21:25])))
colnames(Small) = TABLE.SEVEN.COL.NAMES
colnames(Second) = TABLE.SEVEN.COL.NAMES
colnames(Third) = TABLE.SEVEN.COL.NAMES
colnames(Fourth) = TABLE.SEVEN.COL.NAMES
colnames(Fifth) = TABLE.SEVEN.COL.NAMES


TABLE.SEVEN <- data.frame(rbind(
  Small, Second, Third, Fourth, Fifth
))

colnames(TABLE.SEVEN) = TABLE.SEVEN.COL.NAMES
rownames(TABLE.FIVE) = c("Small", "2", "3", "4", "Big")
print(TABLE.SEVEN)
write.csv(TABLE.FIVE, "T7.CSV")

########################################################### Table 8 ###################################################################
size.rvar.excess <- as.matrix(sweep(data.matrix(size.rvar.filtered[,2:26]), 2, RF, FUN='-'))
MKT.RF = as.matrix(five.factors.filtered[, 2])
SMB = as.matrix(five.factors.filtered[, 3])
HML = as.matrix(five.factors.filtered[, 4])
RMW = as.matrix(five.factors.filtered[, 5])
CMA = as.matrix(five.factors.filtered[, 6])

T8.Panel.A <- matrix(-99, 5, 10)
for(i in 1 : 25) {
  y = as.matrix(size.ni.excess[, i])
  y.x = cbind(y, five.factors.filtered[, 2:6])
  summary = summary(lm(formula = y ~ MKT.RF + SMB + HML, data = y.x))$coefficients
  T8.Panel.A[((i-1) %/% 5 + 1), ((i- 1) %% 5 + 1)] = round(summary[1, 1], 2)
  T8.Panel.A[((i-1) %/% 5 + 1), ((i- 1) %% 5 + 1) + 5] = round(summary[1, 3], 2)
}
rownames(T8.Panel.A) <- c("low", "2", "3", "4", "high")
colnames(T8.Panel.A) <- c("small", "2", "3", "4", "big", "small", "2", "3", "4", "big")
print(T4.Panel.A)

T8.Panel.B <- matrix(-99, 30, 10)
for(i in 1 : 25) {
  y = as.matrix(size.rvar.excess[, i])
  y.x = y.x = cbind(y, five.factors.filtered[, 2:6])
  summary = summary(lm(formula = y ~ MKT.RF + SMB + HML + RMW + CMA, data = y.x))$coefficients
  T8.Panel.B[((i-1) %/% 5 + 1), ((i- 1) %% 5 + 1)] = round(summary[1, 1], 2)
  T8.Panel.B[((i-1) %/% 5 + 1) + 5, ((i- 1) %% 5 + 1)] = round(summary[2, 1], 2)
  T8.Panel.B[((i-1) %/% 5 + 1) + 10, ((i- 1) %% 5 + 1)] = round(summary[3, 1], 2)
  T8.Panel.B[((i-1) %/% 5 + 1) + 15, ((i- 1) %% 5 + 1)] = round(summary[4, 1], 2)
  T8.Panel.B[((i-1) %/% 5 + 1) + 20, ((i- 1) %% 5 + 1)] = round(summary[5, 1], 2)
  T8.Panel.B[((i-1) %/% 5 + 1) + 25, ((i- 1) %% 5 + 1)] = round(summary[6, 1], 2)
  
  T8.Panel.B[((i-1) %/% 5 + 1), ((i- 1) %% 5 + 1) + 5] = round(summary[1, 3], 2)
  T8.Panel.B[((i-1) %/% 5 + 1) + 5, ((i- 1) %% 5 + 1) + 5] = round(summary[2, 3], 2)
  T8.Panel.B[((i-1) %/% 5 + 1) + 10, ((i- 1) %% 5 + 1) + 5] = round(summary[3, 3], 2)
  T8.Panel.B[((i-1) %/% 5 + 1) + 15, ((i- 1) %% 5 + 1) + 5] = round(summary[4, 3], 2)
  T8.Panel.B[((i-1) %/% 5 + 1) + 20, ((i- 1) %% 5 + 1) + 5] = round(summary[5, 3], 2)
  T8.Panel.B[((i-1) %/% 5 + 1) + 25, ((i- 1) %% 5 + 1) + 5] = round(summary[6, 3], 2)
}

colnames(T8.Panel.B) <- c("low", "2", "3", "4", "high", "low", "2", "3", "4", "high")
rownames(T8.Panel.B) <- c("small", "2", "3", "4", "big", "small", "2", "3", "4", "big", 
                          "small", "2", "3", "4", "big", "small", "2", "3", "4", "big", 
                          "small", "2", "3", "4", "big", "small", "2", "3", "4", "big")

print(T8.Panel.B)

# rbind(T8.Panel.A, T8.Panel.B)

TABLE.EIGHT = rbind(
  c("PanelA:Three-factor: Rit–RFt = ai + bi(RMt – RFt) + siSMBt + hiHMLi + eit ", "", "", "", "", "", "", "", "", ""),
  c("", "", "a", "", "", "", "", "t(a)", "", ""), 
  T8.Panel.A[1:5,], 
  c("Panel B: Five-factor: Rit – RFt = ai + bi(RMt – RFt) + siSMBt + hiHMLOt + riRMWt + ciCMAt + eit", "", "", "", "", "", "", "", "", ""),
  c("", "", "a", "", "", "", "", "t(a)", "", ""), 
  T8.Panel.B[1:5,], 
  c("", "", "b", "", "", "", "", "t(b)", "", ""),
  T8.Panel.B[6:10,], 
  c("", "", "s", "", "", "", "", "t(s)", "", ""), 
  T8.Panel.B[11:15,], 
  c("", "", "h", "", "", "", "", "t(h)", "", ""), 
  T8.Panel.B[16:20,], 
  c("", "", "r", "", "", "", "", "t(r)", "", ""), 
  T8.Panel.B[21:25,], 
  c("", "", "c", "", "", "", "", "t(c)", "", ""),
  T8.Panel.B[26:30,]
)
#write.csv(TABLE.EIGHT, "T8.CSV")

########################################################### Table 9 ###################################################################

###########################################
##### 25 Size-AC portfolios #####

TABLE.NINE.COL.NAMES <- c("Low", "2", "3", "4", "High", "Low", "2", "3", "4", "High")

size.ac <- read.csv("25_Portfolios_ME_AC.CSV",header=TRUE)
size.ac <- size.ac[1:column.count,] 
size.ac <- process.input.data(size.ac)
size.ac.filtered <- filter.by.date(size.ac)
size.ac.excess = sweep(size.ac.filtered[,2:26], 2, RF, FUN="-")

size.ac.mean = t(as.matrix(round(apply(size.ac.excess, 2, mean), 2)))
size.ac.sd = t(as.matrix(round(apply(size.ac.excess, 2, sd), 2)))

Small = t(as.matrix(c(size.ac.mean[1, 1:5], size.ac.sd[1, 1:5])))
Second = t(as.matrix(c(size.ac.mean[1, 6:10], size.ac.sd[1, 6:10])))
Third = t(as.matrix(c(size.ac.mean[1, 11:15], size.ac.sd[1, 11:15])))
Fourth = t(as.matrix(c(size.ac.mean[1, 16:20], size.ac.sd[1, 16:20])))
Fifth = t(as.matrix(c(size.ac.mean[1, 21:25], size.ac.sd[1, 21:25])))
colnames(Small) = TABLE.NINE.COL.NAMES
colnames(Second) = TABLE.NINE.COL.NAMES
colnames(Third) = TABLE.NINE.COL.NAMES
colnames(Fourth) = TABLE.NINE.COL.NAMES
colnames(Fifth) = TABLE.NINE.COL.NAMES


TABLE.NINE <- data.frame(rbind(
  Small, Second, Third, Fourth, Fifth
))

colnames(TABLE.NINE) = TABLE.NINE.COL.NAMES
rownames(TABLE.NINE) = c("Small", "2", "3", "4", "Big")
print(TABLE.NINE)
# write.csv(TABLE.NINE, "T9.CSV")

########################################################### Table 10 ###################################################################
size.ac.excess <- as.matrix(sweep(data.matrix(size.ac.filtered[,2:26]), 2, RF, FUN='-'))
MKT.RF = as.matrix(five.factors.filtered[, 2])
SMB = as.matrix(five.factors.filtered[, 3])
HML = as.matrix(five.factors.filtered[, 4])
RMW = as.matrix(five.factors.filtered[, 5])
CMA = as.matrix(five.factors.filtered[, 6])

T10.Panel.A <- matrix(-99, 15, 10)
for(i in 1 : 25) {
  y = as.matrix(size.beta.excess[, i])
  y.x = y.x = cbind(y, five.factors.filtered[, 2:6])
  # regress one
  summary = summary(lm(formula = y ~ MKT.RF + SMB + HML, data = y.x))$coefficients
  T10.Panel.A[((i-1) %/% 5 + 1), ((i- 1) %% 5 + 1)] = round(summary[1, 1], 2)
  T10.Panel.A[((i-1) %/% 5 + 1), ((i- 1) %% 5 + 1) + 5] = round(summary[1, 3], 2)
  
  # regress two
  summary = summary(lm(formula = y ~ MKT.RF + SMB + HML + CMA, data = y.x))$coefficients
  T10.Panel.A[((i-1) %/% 5 + 1) + 5, ((i- 1) %% 5 + 1)] = round(summary[1, 1], 2)
  T10.Panel.A[((i-1) %/% 5 + 1) + 5, ((i- 1) %% 5 + 1) + 5] = round(summary[1, 3], 2)
  
  # regress three
  summary = summary(lm(formula = y ~ MKT.RF + SMB + HML + RMW + CMA, data = y.x))$coefficients
  T10.Panel.A[((i-1) %/% 5 + 1) + 10, ((i- 1) %% 5 + 1)] = round(summary[1, 1], 2)
  T10.Panel.A[((i-1) %/% 5 + 1) + 10, ((i- 1) %% 5 + 1) + 5] = round(summary[1, 3], 2)
}
colnames(T10.Panel.A) <- c("low", "2", "3", "4", "high", "low", "2", "3", "4", "high")
rownames(T10.Panel.A) <- c("small", "2", "3", "4", "big", "small", "2", "3", "4", "big", "small", "2", "3", "4", "big")
print(T10.Panel.A)

T10.Panel.B <- matrix(-99, 30, 10)
for(i in 1 : 25) {
  y = as.matrix(size.beta.excess[, i])
  y.x = y.x = cbind(y, five.factors.filtered[, 2:6])
  # regress one
  summary = summary(lm(formula = y ~ MKT.RF + SMB + HML, data = y.x))$coefficients
  T10.Panel.B[((i-1) %/% 5 + 1), ((i- 1) %% 5 + 1)] = round(summary[4, 1], 2)
  T10.Panel.B[((i-1) %/% 5 + 1), ((i- 1) %% 5 + 1) + 5] = round(summary[4, 3], 2)
  
  # regress two
  summary = summary(lm(formula = y ~ MKT.RF + SMB + HML + CMA, data = y.x))$coefficients
  T10.Panel.B[((i-1) %/% 5 + 1) + 5, ((i- 1) %% 5 + 1)] = round(summary[4, 1], 2)
  T10.Panel.B[((i-1) %/% 5 + 1) + 5, ((i- 1) %% 5 + 1) + 5] = round(summary[4, 3], 2)
  T10.Panel.B[((i-1) %/% 5 + 1) + 10, ((i- 1) %% 5 + 1)] = round(summary[5, 1], 2)
  T10.Panel.B[((i-1) %/% 5 + 1) + 10, ((i- 1) %% 5 + 1) + 5] = round(summary[5, 3], 2)
  
  # regress three
  summary = summary(lm(formula = y ~ MKT.RF + SMB + HML + RMW + CMA, data = y.x))$coefficients
  T10.Panel.B[((i-1) %/% 5 + 1) + 15, ((i- 1) %% 5 + 1)] = round(summary[4, 1], 2)
  T10.Panel.B[((i-1) %/% 5 + 1) + 15, ((i- 1) %% 5 + 1) + 5] = round(summary[4, 3], 2)
  T10.Panel.B[((i-1) %/% 5 + 1) + 20, ((i- 1) %% 5 + 1)] = round(summary[5, 1], 2)
  T10.Panel.B[((i-1) %/% 5 + 1) + 20, ((i- 1) %% 5 + 1) + 5] = round(summary[5, 3], 2)
  T10.Panel.B[((i-1) %/% 5 + 1) + 25, ((i- 1) %% 5 + 1)] = round(summary[6, 1], 2)
  T10.Panel.B[((i-1) %/% 5 + 1) + 25, ((i- 1) %% 5 + 1) + 5] = round(summary[6, 3], 2)
}
colnames(T10.Panel.B) <- c("low", "2", "3", "4", "high", "low", "2", "3", "4", "high")
rownames(T10.Panel.B) <- c("small", "2", "3", "4", "big", "small", "2", "3", "4", "big", 
                          "small", "2", "3", "4", "big", "small", "2", "3", "4", "big", 
                          "small", "2", "3", "4", "big", "small", "2", "3", "4", "big")

print(T10.Panel.B)

TABLE.TEN = rbind(
  c("PanelA: Regression intercepts", "", "", "", "", "", "", "", "", ""),
  c("", "", "a", "", "", "", "", "t(a)", "", ""), 
  c("Three-factor: Mkt, SMB, and HML", "", "", "", "", "", "", "", "", ""),
  T10.Panel.A[1:5,], 
  c("Four-factor: Mkt, SMB, HML, and CMA", "", "", "", "", "", "", "", "", ""),
  T10.Panel.A[6:10,], 
  c("Five-factor: Mkt, SMB, HMLO, RMW, and CMA", "", "", "", "", "", "", "", "", ""),
  T10.Panel.A[11:15,], 
  c("Panel B: Regression slopes", "", "", "", "", "", "", "", "", ""),
  c("Three-factor: Rit – RFt = ai + bi(RMt – RFt) + siSMBt + hiHMLt + eit", "", "", "", "", "", "", "", "", ""),
  c("", "", "h", "", "", "", "", "t(h)", "", ""), 
  T10.Panel.B[1:5,], 
  c("Four-factor: Rit – RFt = ai + bi(RMt – RFt) + siSMBt + hiHMLt + ciCMAt + eit", "", "", "", "", "", "", "", "", ""),
  c("", "", "h", "", "", "", "", "t(h)", "", ""),
  T10.Panel.B[6:10,], 
  c("", "", "c", "", "", "", "", "t(c)", "", ""), 
  c("Five-factor: Rit – RFt = ai + bi(RMt – RFt) + siSMBt + hiHMLOt + riRMWt + ciCMAt + eit", "", "", "", "", "", "", "", "", ""),
  T10.Panel.B[11:15,], 
  c("", "", "h", "", "", "", "", "t(h)", "", ""), 
  T10.Panel.B[16:20,], 
  c("", "", "r", "", "", "", "", "t(r)", "", ""), 
  T10.Panel.B[21:25,], 
  c("", "", "c", "", "", "", "", "t(c)", "", ""),
  T10.Panel.B[26:30,]
)
#write.csv(TABLE.TEN, "T10.CSV")