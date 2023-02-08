#source("helpers/kdepairs.default.R")
#source("helpers/kdepairs.R")
source("helpers/kde_ggpairs.R")

library(tools)
library(bayestestR)

cat("Did you delete all the variables you wanted? If no, run 'rm(list=ls())'")
#rm(list=ls())
# ----------------------------------------------------------------------
# name of the data file

prf<- "../TREATMENT_CALIBRATIONS_RUN_3/CAL_KSC_PAC_AE_EXP_IC_SET_2_RUN_4_NO_R/curgen_db_0"
it <- c("00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19")

for(i in 1:20)
{
#i=14
fname <- paste0(prf,it[i],".txt")  


# names of variables
l = c(expression(D[u]),expression(s),expression(k),expression(chi),expression(D[f]),expression(sigma ^2)) #,expression(d),expression(chi),expression(D[f]),expression(r)


discrepancy <- 0  # is the last column of the data discrepancy of likelihood?
truelik <- 1  # use likelihood or log-likelihood for coloring

# ----------------------------------------------------------------------
print_tab <- function(text, x, nd, w) {
  cat(text, formatC(x, digits=nd, format='e', width=w), "\n", sep="\t")
}

print_nd <- function(x, nd) {
  cat(format(round(x, digits=nd), nsmall=nd))
}

cat("Reading file", fname,"\n")
data <- read.table(fname)
data <- array(data=unlist(data), dim=dim(data))
nd <- dim(data)[2]  # dimension of the samples

# mean
md <- colMeans(data[, 1:nd-1])
print_tab("means:\t", md, 4, 8)

# standard deviation
m2d <- colMeans(data[, 1:nd-1]^2)
sd <- sqrt(m2d-md^2)
print_tab("stds:\t", sd, 4, 8)

# most probable parameters
best_id <- which.max(data[, nd-1])
best <- data[best_id, 1:nd-1]
print_tab("MPVs:\t", best, 4, 8)

# quantiles
q1 <- c(); for (i in seq(1, nd-1)) q1[i] <- quantile(data[, i], probs=c(0.05))
print_tab("q-0.05:  ", q1, 4, 8)
q2 <- c(); for (i in seq(1, nd-1)) q2[i] <- quantile(data[, i], probs=c(0.95))
print_tab("q-0.95:  ", q2, 4, 8)


# print table
gname <- paste0(file_path_sans_ext(fname), "_tmcmc_params.txt")
sink(gname)
print_tab("means:", md, 4, 8)
print_tab("stds:", sd, 4, 8)
print_tab("MPVs:", best, 4, 8)
print_tab("q+0.05:", q1, 4, 8)
print_tab("q-0.95:", q2, 4, 8)
sink()

bname <- paste0(file_path_sans_ext(fname), "_best.txt")
#sink(bname)
#best
#sink()
for(i in 1:6){
  
  cat(best[i],file=bname,sep="\n",append = TRUE)
}

#write(best,bname)

dat <- data.frame(data)
colnames(dat) = l
ci_eti <- ci(dat,ci=0.95, method = "HDI")
fname1 <-paste0(file_path_sans_ext(fname), "_CI95.txt")
sink(fname1)
for(i in (3:nd-2))
  print_tab(ci_eti$Parameter[i], c(ci_eti$CI_low[i], ci_eti$CI_high[i]), 4,8)
sink()

# print for latex table
a <- c(); digits <- 3
for (i in seq(1, nd-1)) {
    t <- capture.output(cat(print_nd(md[i], digits)))
    a <- paste(c(a, t), collapse = "")
    a <- paste(c(a, " & ["), collapse = "")
    t <- capture.output(cat(print_nd(q1[i], digits)))
    a <- paste(c(a, t), collapse = "")
    a <- paste(c(a, ", "), collapse = "")
    t <- capture.output(cat(print_nd(q2[i], digits)))
    a <- paste(c(a, t), collapse = "")
    a <- paste(c(a, "] & "), collapse = "")
}

#write.table(a, file=fres,sep = "\t", row.names = T)
a <- capture.output(cat(substr(a, 1, nchar(a)-3)))
a <- paste("latex:", a)
a <- paste(a, "\\\\ \n")
cat(a)



if(discrepancy) data[, nd] <-    -data[, nd]
if(truelik)     data[, nd] <- exp(data[, nd])

pname <- paste0(file_path_sans_ext(fname), ".png")
png(pname, width=2500, height=2500, units='px', res=300, pointsize=15, type="cairo", family="times")
pm<-kdeggpairs(dat, n_1d=20, n_2d=200, labels=l) #n_2d=200
print(pm)
#kdepairs(data, n_1d=20, n_2d=200, labels=l) #n_2d=200
dev.off()


}
