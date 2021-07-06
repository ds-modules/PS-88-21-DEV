#Replication Code for Taylor C. Boas and F. Daniel Hidalgo, "Controlling the Airwaves: Incumbency Advantage and Community Radio in Brazil," American Journal of Political Science

#NOTE: Comments indicate where the code begins to recreate each table (or plot) in the main text and in the supporting information. The code below produces latex code for each table, but note that we reformatted the tables in some instances so they may appear somewhat different in the paper. 

#Analysis conducted in R version 2.12.2

#install the following packages
library(car)
library(xtable)
library(Matching)
library(reporttools)

#Set working directory in the following line:
setwd("")
load("replication.RData")

####
#Create Helper Functions For Estimation and Plots
###
###Discontinuity Plots###
reg.discont.plot <- function(x,force.var,treat,bins=50, xlab="Forcing Variable",ylab="Dependent Variable", bw = .5,SE=FALSE,y.max=.9, y.min=.1,main=""){
	library(Hmisc)
	
	not.missing <- which(is.na(x)==FALSE & is.na(force.var)==FALSE & is.na(treat)==FALSE)
	x <- x[not.missing]
	force.var <- force.var[not.missing]
	treat <- treat[not.missing]

	bins.controls <- levels(cut2(force.var[treat==0],g=round(bins)))

	bin.midpoints.controls <- (as.numeric(as.numeric(sub(pattern= "\\[(.*),.*",x= bins.controls,replacement= "\\1"))) +  as.numeric(sub(pattern= ".*,(.*).",x= bins.controls,replacement="\\1")))/2
	
	bins.treated <- levels(cut2(force.var[treat==1],g=round(bins)))
	bin.midpoints.treated <- (as.numeric(as.numeric(sub(pattern= "\\[(.*),.*",x= bins.treated,replacement= "\\1"))) +  as.numeric(sub(pattern= ".*,(.*).",x= bins.treated,replacement="\\1")))/2

	#mean within bins for controls, use the "cut " function
	bin.means.controls <- tapply(x[treat==0],cut2(force.var[treat==0],g=round(bins)),mean,na.rm=TRUE)
	bin.means.treated <- tapply(x[treat==1],cut2(force.var[treat==1],g=round(bins)),mean,na.rm=TRUE)
   
	y.min <- quantile(c(bin.means.controls,bin.means.treated),probs=y.min,na.rm=TRUE)
	y.max <- quantile(c(bin.means.controls,bin.means.treated),probs=y.max,na.rm=TRUE)
	x.min <- min(force.var,na.rm=TRUE)
	x.max <- max(force.var,na.rm=TRUE)

	plot(1,type="n",ylim=c(y.min,y.max),xlim=c(x.min,x.max),xlab=xlab, ylab=ylab,main=main)


	points(bin.midpoints.controls, bin.means.controls,pch=19,col="blue")
	points(bin.midpoints.treated,bin.means.treated,pch=19,col="red")
	abline(v=0,lwd=1,lty=3)
	

	local.poly.control <- loess(x[treat==0]~force.var[treat==0],span=bw,degree=1)
	predict.local.poly.control <- predict(local.poly.control,sort(force.var[treat==0]),se=TRUE)
	lines(sort(force.var[treat==0]),predict.local.poly.control$fit, col="blue",lwd=2)
	if(SE==TRUE){
	lines(sort(force.var[treat==0]),predict.local.poly.control$fit+(1.96*predict.local.poly.control$se.fit), col="blue",lwd=2,lty=2)
	lines(sort(force.var[treat==0]),predict.local.poly.control$fit-(1.96*predict.local.poly.control$se.fit), col="blue",lwd=2,lty=2)
}
	
	local.poly.treat <- loess(x[treat==1]~force.var[treat==1],span=bw,degree=1,)
	predict.local.poly.treat <- predict(local.poly.treat,sort(force.var[treat==1]),se=TRUE)
	lines(sort(force.var[treat==1]),predict.local.poly.treat$fit, col="red",lwd=2)
	if(SE==TRUE){
	lines(sort(force.var[treat==1]),predict.local.poly.treat$fit+(1.96*predict.local.poly.treat$se.fit), col="red",lwd=2,lty=2)
	lines(sort(force.var[treat==1]),predict.local.poly.treat$fit-(1.96*predict.local.poly.treat$se.fit), col="red",lwd=2,lty=2)
}
	

}

##
###Function for Bandwidth Selection with Cross-Validation###
##
cross.valid.bw <- function(y,force.var,force.var.window.left,force.var.window.right,cut.point,h.length=50,min.h=.05,max.h=.5){
  not.missing <- is.na(y)==FALSE & is.na(force.var)==FALSE
  y <- y[not.missing]
  force.var <- force.var[not.missing]

  h <- seq(from = quantile(abs(force.var),probs=min.h), to = quantile(abs(force.var),probs=max.h), length.out=h.length)


  y.window <- y[force.var>force.var.window.left & force.var<force.var.window.right]
  force.var.window <- force.var[force.var>force.var.window.left & force.var<force.var.window.right]

  cv <- rep(NA,times=h.length)
  for(j in 1:length(h)){
    y.hat <- rep(NA,times=length(y.window))
    for (i in 1:length(y.window)){
      c <- force.var.window[i]
      if(c<cut.point){
        x.bw <- force.var[force.var< c & force.var>(c- h[j])]
        y.bw <- y[force.var< c & force.var>(c- h[j])]
      }
      if(c>=cut.point){
        x.bw <- force.var[force.var> c & force.var<(c+ h[j])]
        y.bw <- y[force.var> c & force.var<(c+ h[j])]
      }
      if(length(x.bw)==0){
        next
      }
      x.bw <- x.bw 
      y.hat[i] <- coef(lm(y.bw~I(x.bw-c)))[1]
    }
    cv[j] <- sum((y.window - y.hat)^2,na.rm=TRUE)/length(y.window[is.na(y.window)==FALSE])
    print(j)
    print(cv[j])
  }
  h.opt <- h[which(cv==min(cv))]
  return(h.opt)
}

#Function for Discontinuity Estimates
discont.point.est <- function(y, force.var, treat, bw.loc.lin, y.name="dv", bw.discont.samp = c(40, 20, 10)){
  library(car)
  local.lin.cv <- lm(y~treat*force.var, subset= abs(force.var) <= bw.loc.lin)
  discont.sample.1 <- lm(y~treat,subset=abs(force.var)<=bw.discont.samp[1])
  discont.sample.2 <- lm(y~treat,subset=abs(force.var)<=bw.discont.samp[2])
  discont.sample.3 <- lm(y~treat,subset=abs(force.var)<=bw.discont.samp[3])
  ols <- lm(y~treat)
  specifications <- c("Local Linear Regression", paste("Discontinuity Sample (vote margin < ", bw.discont.samp[1],")", sep =""), paste("Discontinuity Sample (vote margin < ", bw.discont.samp[2],")", sep =""), paste("Discontinuity Sample (vote margin < ", bw.discont.samp[3],")", sep =""),"Full Sample")
  coefs <- c(coef(local.lin.cv)[2],   
coef(discont.sample.1)[2], coef(discont.sample.2)[2],  coef(discont.sample.3)[2],  coef(ols)[2])  
  se <- c(sqrt(hccm(local.lin.cv)[2,2]), sqrt(hccm(discont.sample.3)[2,2]), sqrt(hccm(discont.sample.2)[2,2]), sqrt(hccm(discont.sample.1)[2,2]), sqrt(hccm(ols)[2,2]))
  n <- c(length(y[abs(force.var)<=bw.loc.lin]), length(y[abs(force.var)<=bw.discont.samp[1]]), length(y[abs(force.var)<=bw.discont.samp[2]]), length(y[abs(force.var)<=bw.discont.samp[3]]),  length(y))
  results <- data.frame(variable = as.factor(rep(y.name,times=5)),specification = as.factor(specifications), coefs = coefs, se = se,n=n)
  results$specification <- factor(results$specification, c("Full Sample","Local Linear Regression", paste("Discontinuity Sample (vote margin < ", bw.discont.samp[1],")", sep =""), paste("Discontinuity Sample (vote margin < ", bw.discont.samp[2],")", sep =""), paste("Discontinuity Sample (vote margin < ", bw.discont.samp[3],")", sep ="")))
  return(results)
}

#####
#TABLES AND FIGURES FROM PAPER
#####
##
#Table 1
#
h.cv <- cross.valid.bw(y = applied.data$post.approved, force.var  = applied.data$raw.vote.margin,force.var.window.right=max(abs(applied.data$raw.vote.margin)),force.var.window.left=-1*max(abs(applied.data$raw.vote.margin)),cut.point=0,h.length=100,min.h=.01,max.h=.6)
#h.cv <-  165.1333


#Check Balance
#Decision to apply
apply.check <- with(full.data, discont.point.est(y = post.applied, force.var = raw.vote.margin, treat = elected, bw.loc.lin = h.cv, y.name = "Apply"))
#Prior approval
prior.approval.check <-  with(full.data, discont.point.est(y = pre.approved, force.var = raw.vote.margin, treat = elected, bw.loc.lin = h.cv, y.name = "Prior Approval"))
#other covariates
yob.check <- with(applied.data,discont.point.est(yob,raw.vote.margin,elected,h.cv,"Year of Birth"))
electorate.check <- with(applied.data,discont.point.est(log(total.votes),raw.vote.margin,elected,h.cv,"Log Electorate"))
coalition.check <- with(applied.data,discont.point.est(log(coal.votes),raw.vote.margin,elected,h.cv,"Log Coalition Votes"))
edu.primary.check <- with(applied.data,discont.point.est(edu_primary,raw.vote.margin,elected,h.cv,"Primary Education"))
pre.approved.check <- with(applied.data,discont.point.est(pre.approved,raw.vote.margin,elected,h.cv,"Approved Before Election"))
uf.mg.check <- with(applied.data,discont.point.est(uf_mg,raw.vote.margin,elected,h.cv,"Minas Gerais"))
uf.sp.check <- with(applied.data,discont.point.est(uf_sp,raw.vote.margin,elected,h.cv,"Sao Paulo"))
occ.agricultor.check  <-with(applied.data,discont.point.est(occ_agricultor,raw.vote.margin,elected,h.cv,"Occupation: Agriculture"))
occ.business.check  <-with(applied.data,discont.point.est(occ_business,raw.vote.margin,elected,h.cv,"Occupation: Business"))
pmdb.check  <-with(applied.data,discont.point.est(party_pmdb,raw.vote.margin,elected,h.cv,"Party: PMDB"))
psdb.check  <-with(applied.data,discont.point.est(party_psdb,raw.vote.margin,elected,h.cv,"Party: PSDB"))
pfl.check  <-with(applied.data,discont.point.est(party_pfl,raw.vote.margin,elected,h.cv,"Party: PFL"))
pt.pres.98.check <-with(applied.data,discont.point.est(pt_pres_1998,raw.vote.margin,elected,h.cv,"Lula Vote Share, 1998"))
gini.check  <-with(applied.data,discont.point.est(gini_2000,raw.vote.margin,elected,h.cv,"Gini (2000)"))
latitude.check  <-with(applied.data,discont.point.est(latitude,raw.vote.margin,elected,h.cv,"Latitude"))
longitude.check  <-with(applied.data,discont.point.est(longitude,raw.vote.margin,elected,h.cv,"Longitude"))
hdi.check  <-with(applied.data,discont.point.est(hdi_2000,raw.vote.margin,elected,h.cv,"HDI (2000)"))
income.check  <-with(applied.data,discont.point.est(income_2000,raw.vote.margin,elected,h.cv,"GDP per Cap (2000)"))
pt.pref.check  <-with(applied.data,discont.point.est(pt_pref_2000,raw.vote.margin,elected,h.cv,"PT Mayor (2000)"))
time.since.app.check <- with(applied.data,discont.point.est(time.since.app,raw.vote.margin,elected,h.cv,"Time Since Application"))

#Put together into one data.set and produce a latex table
bal.check.results <- data.frame(rbind(apply.check, prior.approval.check,yob.check,coalition.check,edu.primary.check,uf.sp.check,uf.mg.check,occ.agricultor.check,occ.business.check,pmdb.check,psdb.check,pfl.check,pt.pres.98.check,gini.check,latitude.check,longitude.check,hdi.check,pt.pref.check,electorate.check,time.since.app.check))
n.var <- length(unique(bal.check.results$variable))
bal.check.table <- (matrix(nrow=n.var*2,ncol=6))
bal.check.table[seq(1,n.var*2,by=2),1] <- unique(as.character(bal.check.results$variable))
bal.check.table[seq(2,n.var*2,by=2),1] <- ""
bal.check.table[seq(1,n.var*2,by=2),2] <- bal.check.results$coef[bal.check.results$specification=="Full Sample"]
bal.check.table[seq(1,n.var*2,by=2),3] <- bal.check.results$coef[bal.check.results$specification=="Local Linear Regression"]
bal.check.table[seq(1,n.var*2,by=2),4] <- bal.check.results$coef[bal.check.results$specification=="Discontinuity Sample (vote margin < 40)"]
bal.check.table[seq(1,n.var*2,by=2),5] <- bal.check.results$coef[bal.check.results$specification=="Discontinuity Sample (vote margin < 20)"]
bal.check.table[seq(1,n.var*2,by=2),6] <- bal.check.results$coef[bal.check.results$specification=="Discontinuity Sample (vote margin < 10)"]
bal.check.table[seq(2,n.var*2,by=2),2] <- bal.check.results$se[bal.check.results$specification=="Full Sample"]
bal.check.table[seq(2,n.var*2,by=2),3] <- bal.check.results$se[bal.check.results$specification=="Local Linear Regression"]
bal.check.table[seq(2,n.var*2,by=2),4] <- bal.check.results$se[bal.check.results$specification=="Discontinuity Sample (vote margin < 10)"]
bal.check.table[seq(2,n.var*2,by=2),5] <- bal.check.results$se[bal.check.results$specification=="Discontinuity Sample (vote margin < 20)"]
bal.check.table[seq(2,n.var*2,by=2),6] <- bal.check.results$se[bal.check.results$specification=="Discontinuity Sample (vote margin < 40)"]
colnames(bal.check.table) <- c("Variable","Full Sample", "Loc. Lin. Reg", "Discont Sample (x<40)","Discont Sample (x<20)","Discont Sample (x<10)")
bal.check.table[,2:6] <- apply(bal.check.table[,2:6],2,function(x) round(as.numeric(x),digits=2))
bal.check.table[seq(2,n.var*2,by=2),2:6] <- paste("(",bal.check.table[seq(2,n.var*2,by=2),2:6],")",sep="")
bal.check.table <- rbind(bal.check.table[1:4,], c("n",unique(bal.check.results$n[bal.check.results$specification=="Full Sample"][1:2]),unique(bal.check.results$n[bal.check.results$specification=="Local Linear Regression"][1:2]),unique(bal.check.results$n[bal.check.results$specification=="Discontinuity Sample (vote margin < 40)"][1:2]),unique(bal.check.results$n[bal.check.results$specification=="Discontinuity Sample (vote margin < 20)"][1:2]),unique(bal.check.results$n[bal.check.results$specification=="Discontinuity Sample (vote margin < 10)"][1:2])), bal.check.table[5:nrow(bal.check.table),])
bal.check.table[1,2:6] <- as.character(signif(apply.check[c(5,1:4),3],2))
bal.check.table[2,2:6] <- paste("(",as.character(signif(apply.check[c(5,1:4),4]),2),")",sep="")
bal.check.table[3,2:6] <- as.character(signif(prior.approval.check[c(5,1:4),3],2))
bal.check.table[4,2:6] <- paste("(",as.character(signif(prior.approval.check[c(5,1:4),4]),2),")",sep="")

bal.check.table <- xtable(bal.check.table,display=c("s","s",rep("fg",5)), caption = "Balance Statistics")
digits(bal.check.table) <-2
align(bal.check.table) <- "rr|rrrrr"
print(bal.check.table,latex.environment = "center", include.rownames = FALSE,hline.after=c(0,nrow(bal.check.table)-1),floating.environment = "sidewaystable")

##
#TABLE 2
##
loc.lin.app <- (lm(post.approved~elected*raw.vote.margin,data=applied.data,subset=abs(raw.vote.margin)<=h.cv))
loc.lin.rej<- (lm(1-pending_yes~elected*raw.vote.margin,data=applied.data,subset=abs(raw.vote.margin)<=h.cv))
discont40.app <- lm(post.approved~elected, data=applied.data,subset = abs(raw.vote.margin)<=40)
discont40.rej <- lm(1-pending_yes~elected, data=applied.data,subset = abs(raw.vote.margin)<=40)
discont20.app <- lm(post.approved~elected, data=applied.data,subset = abs(raw.vote.margin)<=20)
discont20.rej <- lm(1-pending_yes~elected, data=applied.data,subset = abs(raw.vote.margin)<=20)
discont10.app <- lm(post.approved~elected, data=applied.data,subset = abs(raw.vote.margin)<=10)
discont10.rej <- lm(1-pending_yes~elected, data=applied.data,subset = abs(raw.vote.margin)<=10)

app.coefs <- c(coef(loc.lin.app)[2],coef(discont40.app)[2],coef(discont20.app)[2],coef(discont10.app)[2])
app.se <- sqrt(c(vcov(loc.lin.app)[2,2],vcov(discont40.app)[2,2],vcov(discont20.app)[2,2],vcov(discont10.app)[2,2]))
app.n <- c(nrow(loc.lin.app$model),nrow(discont40.app$model),nrow(discont20.app$model), nrow(discont10.app$model))
app.table <- (matrix(rbind(app.coefs,app.se,app.n),nrow=3))

rej.coefs <- c(coef(loc.lin.rej)[2],coef(discont40.rej)[2],coef(discont20.rej)[2],coef(discont10.rej)[2])
rej.se <- sqrt(c(vcov(loc.lin.rej)[2,2],vcov(discont40.rej)[2,2],vcov(discont20.rej)[2,2],vcov(discont10.rej)[2,2]))
rej.n <- c(nrow(loc.lin.rej$model),nrow(discont40.rej$model),nrow(discont20.rej$model), nrow(discont10.rej$model))
rej.table <- (matrix(rbind(rej.coefs,rej.se,rej.n),nrow=3))

xtable(rbind(app.table,rej.table))

##
#Table 3
##
match.pctVV <- Match(Y= match.data$pctVV, Tr = match.data$treat, X = genmatch.covar, Weight.matrix=genmatch.output, BiasAdjust=FALSE, caliper = c(rep(10,43),.5))
match.elected <- Match(Y= match.data$elected, Tr = match.data$treat, X = genmatch.covar, Weight.matrix=genmatch.output,  caliper = c(rep(10,43),.5), BiasAdjust=FALSE)
bal.data <- with(match.data, data.frame(match.party, male, occ, edu, yob, lat, long, ran.prior, incumbent, log.valid.votes, prior.pctVV, party.prior.pctVV, elec.year, log.total.assets, gini_2000, pt_pres_1998, income_2000, psdb_2000, pt_2000, hdi_2000, uf.ba, uf.sp, uf.mg, uf.rs, log.num.apps))
bal.fmla <- as.formula(paste("treat~",paste(names(bal.data),collapse="+")))
match.bal <- MatchBalance(bal.fmla, data = data.frame(treat = match.data$treat, bal.data), match.out = match.pctVV,nboots=3000)
before.sd.diff <- sapply(match.bal$BeforeMatching, function(x) x[[1]])
after.sd.diff <- sapply(match.bal$AfterMatching, function(x) x[[1]])
before.t.p <- matrix(sapply(match.bal$BeforeMatching, function(x) x$tt[3]))
after.t.p <- matrix(sapply(match.bal$AfterMatching, function(x) x$tt[3]))
before.ks.p <- matrix(sapply(match.bal$BeforeMatching, function(x) x$ks$ks.boot.pvalue))
after.ks.p <- matrix(sapply(match.bal$AfterMatching, function(x) x$ks$ks.boot.pvalue))
bal.stats <- data.frame(var = names(data.frame(model.matrix(as.formula(paste("~",paste(names(bal.data),collapse="+"))), bal.data)[,-1])), before.sd.diff, after.sd.diff, before.t.p = unlist(before.t.p), after.t.p = unlist(after.t.p), before.ks.p, after.ks.p)
bal.stats <- bal.stats[c(4,6,12,15,18:27,29:49),]
rownames(bal.stats)  <- c("Party: PFL", "Party: PMDB", "Party: PSDB", "Party: PT", "Male", "Occupation: Blue Collar", "Occupation: Education", "Occupation: Government", "Occupation: Media", "Occupation: None", "Occupation: Other", "Occupation: Politician", "Occupation: White Collar", "Education: Some Superior or More", "Year of Birth", "Latitude", "Longitude", "Ran Previously", "Incumbency", "Log Electorate", "Prior Vote Share", "Party Prior Vote Share", "Election Year", "Total Assets", "2000 Gini", "PT Pres Vote Share (1998)", "GDP Per Capita (2000)", "PSDB Mayor Vote Share (2000)", "PT Mayor Vote Share (2000)", "HDI (2000)", "State: Bahia", "State: Sao Paulo", "State: Minas Gerais", "State: Rio Grande do Sul", "Log Number of Applications")
bal.stats$before.ks.p[bal.stats$before.ks.p=="NULL"] <- NA
bal.stats$after.ks.p[bal.stats$after.ks.p=="NULL"] <- NA
bal.stats$before.ks.p <- unlist(bal.stats$before.ks.p)
bal.stats$after.ks.p <- unlist(bal.stats$after.ks.p)
bal.stats <- bal.stats[order(abs(bal.stats$before.sd.diff), decreasing = TRUE),]
bal.stats <- bal.stats[,-1]
bal.check.table <- xtable(bal.stats, display=c("s",rep("fg",6)), caption = "Balance Statistics")
digits(bal.check.table) <-2
align(bal.check.table) <- "r|rrrrrr"
print(bal.check.table,latex.environment = "center", include.rownames = TRUE,hline.after=c(0,nrow(bal.check.table)),floating.environment = "sidewaystable")

##
#Table 4
##

results.tab <- data.frame("Pct Valid Votes" = c((match.pctVV$est), (match.pctVV$se), (length(match.pctVV$mdata$Y))),
"Elected" = c((match.elected$est), (match.elected$se), (length(match.elected$mdata$Y))))
rownames(results.tab) <- c( "ATT Estimate", "SE", "n")
results.tab <- xtable( results.tab, display=c("s",rep("fg",2)), caption = "Estimates")
digits(results.tab) <-2
align(results.tab) <- "r|rr"
print(results.tab,latex.environment = "center", include.rownames = TRUE,hline.after=c(0))
##
#Figure 1
##
bins <- 15
pdf("approved_rejected_discont.pdf", height = 11, width  = 5)#use this command to create pdf
par(mfrow=c(2,1))
reg.discont.plot(applied.data$post.approved[abs(applied.data$raw.vote.margin)<250], applied.data$raw.vote.margin[abs(applied.data$raw.vote.margin)<250], treat = applied.data$elected[abs(applied.data$raw.vote.margin)<250],bins=15, bw=1, SE=TRUE, y.max=1, y.min=0, xlab=expression(paste("Vote Margin, ",M[ij])), ylab= "Application Approved After Election", main= "Application Approved After Election")
reg.discont.plot(1-applied.data$pending_yes[abs(applied.data$raw.vote.margin)<250], applied.data$raw.vote.margin[abs(applied.data$raw.vote.margin)<250], treat=applied.data$elected[abs(applied.data$raw.vote.margin)<250],bins=15,bw=1, SE=TRUE,y.max=1,y.min=0, xlab=expression(paste("Vote Margin, ",M[ij])), ylab="Application Rejected After Election", main="Application Rejected After Election")
dev.off()


###
#TABLES AND FIGURES FROM SUPPLEMENTARY APPENDIX
###
##
#Table 1
##
discont.desc.vars <- c("abs.vote.margin", "coal.votes", "total.votes", "gini_2000", "hdi_2000", "psdb_2000", "pt_pres_1998", "time.since.app", "yob", "edu_primary", "elected", "party_pmdb", "party_pt", "uf_mg", "occ_bureaucrat", "occ_business", "uf_sp", "party_psdb", "party_pfl", "occ_agricultor", "uf_BA", "pt_pref_2000")
discont.desc.cap <- "Effect of Incumbency: Descriptive Statistics"

applied.data$abs.vote.margin <- abs(applied.data$raw.vote.margin)


discont.desc <- applied.data[,discont.desc.vars]
names(discont.desc) <- c("Abs(Vote Margin)", "Coalition Votes", "Electorate", "Income Gini (2000)", "HDI (2000)", "PSDB (2000)", "PT President Vote Share (1998)", "Time Since Application", "Year of Birth","Primary Education", "Elected", "Party: PMDB", "Party: PT", "State: MG", "Occupation: Bureaucrat", "Occupation: Business", "State: SP", "Party: PSDB", "Party: PFL", "Occupation: Agriculture", "State: BA", "PT Mayor (2000)" )
discont.desc <- discont.desc[,order(names(discont.desc))]
discont.desc.table <-tableContinuous(vars = discont.desc, cap = discont.desc.cap, lab = "tab:discont_desc", longtable=FALSE, stats = c("n", "min", "median", "mean", "max", "s"), prec=2)
##
#Table 2
##
match.data$occ.media <- ifelse(match.data$occ=="Media",1,0)
match.data$occ.bluecollar <- ifelse(match.data$occ=="Blue collar",1,0)
match.data$occ.whitecollar <- ifelse(match.data$occ=="White collar",1,0)
match.data$occ.none <- ifelse(match.data$occ=="None",1,0)
match.data$occ.politician <- ifelse(match.data$occ=="Politician",1,0)
match.data$occ.government <- ifelse(match.data$occ=="Government",1,0)
match.data$occ.other <- ifelse(match.data$occ=="Other",1,0)
match.data$occ.education <- ifelse(match.data$occ=="Education",1,0)
match.data$edu.some.superior <- ifelse(match.data$edu == "Some Superior or More", 1, 0)
match.data$party.pt <- ifelse(match.data$party=="PT", 1, 0)
match.data$party.psdb <- ifelse(match.data$party=="PSDB", 1, 0)
match.data$party.pmdb <- ifelse(match.data$party=="PMDB", 1, 0)
match.data$party.pfl <- ifelse(match.data$party=="PFL", 1, 0)
match.data$total.assets <- exp(match.data$log.total.assets)-1


match.desc.vars <- c("uf.rs", "occ.media", "male","occ.bluecollar", "occ.whitecollar", "uf.mg", "edu.some.superior","uf.ba", "occ.none", "party.pmdb", "incumbent","party.psdb","party.pfl","uf.sp","occ.politician", "occ.government","ran.prior","occ.other","occ.education", "num.apps","pt_pres_1998", "total.assets","valid.votes","pt_2000", "income_2000", "hdi_2000", "yob", "lat", "party.prior.pctVV", "psdb_2000", "prior.pctVV","long")
match.desc.cap <- "Effect of Control: Descriptive Statistics"
match.desc <- match.data[,match.desc.vars]
names(match.desc) <- c("State: Rio Grande do Sul", "Occupation: Media", "Male", "Occupation: Blue Collar", "Occupation: White Collar", "State: Minas Gerais", "Education: Some Superior or More", "State: Bahia", "Occupation: None", "Party:PMDB", "Incumbent", "Party: PSDB", "Party: PFL", "UF: SP", "Occupation: Politician", "Occupation: Government", "Ran Previously", "Occupation: Other", "Occupation: Education", "Number of Applications", "PT Presidential Vote Share (1998)", "Total Asset Value", "Electorate", "PT Mayoral Vote Share (2000)", "GDP per capita (2000)", "HDI (2000)", "Year of Birth", "Latitude", "Party's Prior Vote Share", "PSDB Mayoral Vote Share (2000)", "Prior Vote Share", "Longitude")
match.desc <- match.desc[,order(names(match.desc))]
tableContinuous(vars = match.desc, cap = match.desc.cap, lab = "tab:match_desc", longtable=FALSE, stats = c("n", "min", "median", "mean", "max", "s"), prec=2)

##
#Table 3
##
#Get cross-validation bandwidth
h.cv.elec <- cross.valid.bw(y = applied.data$post.approved, force.var  = applied.data$vote.margin.electorate,force.var.window.righ=max(abs(applied.data$vote.margin.electorate)),force.var.window.left=-1*max(abs(applied.data$vote.margin.electorate)),cut.point=0,h.length=100,min.h=.01,max.h=.6)
#h.cv.elec <- 0.009737304

approved.results.elec <- with(applied.data,discont.point.est(y = post.approved, force.var = vote.margin.electorate,treat = elected,bw.loc.lin= h.cv.elec, y.name ="License Approval", bw.discont.samp = c(.006, .004, .002)))
rejected.results.elec <- with(applied.data,discont.point.est(y = 1-pending_yes, force.var = vote.margin.electorate,treat = elected,bw.loc.lin= h.cv.elec, y.name ="License Rejection", bw.discont.samp = c(.006, .004, .002)))

##Covariates##
#Decision to apply
apply.check.elec <- with(full.data, discont.point.est(y = post.applied, force.var = vote.margin.electorate, treat = elected, bw.loc.lin = h.cv.elec, y.name = "Apply", bw.discont.samp = c(.006, .004, .002)))
#Prior approval
prior.approval.check.elec <-  with(full.data, discont.point.est(y = pre.approved, force.var = vote.margin.electorate, treat = elected, bw.loc.lin = h.cv.elec, y.name = "Prior Approval", bw.discont.samp = c(.006, .004, .002)))
#other covariates
yob.check.elec <- with(applied.data,discont.point.est(yob,vote.margin.electorate,elected,h.cv.elec,"Year of Birth", bw.discont.samp = c(.006, .004, .002)))
electorate.check.elec <- with(applied.data,discont.point.est(log(total.votes),vote.margin.electorate,elected,h.cv.elec,"Log Electorate", bw.discont.samp = c(.006, .004, .002)))
coalition.check.elec <- with(applied.data,discont.point.est(log(coal.votes),vote.margin.electorate,elected,h.cv.elec,"Log Coalition Votes", bw.discont.samp = c(.006, .004, .002)))
edu.primary.check.elec <- with(applied.data,discont.point.est(edu_primary,vote.margin.electorate,elected,h.cv.elec,"Primary Education", bw.discont.samp = c(.006, .004, .002)))
pre.approved.check.elec <- with(applied.data,discont.point.est(pre.approved,vote.margin.electorate,elected,h.cv.elec,"Approved Before Election", bw.discont.samp = c(.006, .004, .002)))
uf.mg.check.elec <- with(applied.data,discont.point.est(uf_mg,vote.margin.electorate,elected,h.cv.elec,"Minas Gerais", bw.discont.samp = c(.006, .004, .002)))
uf.sp.check.elec <- with(applied.data,discont.point.est(uf_sp,vote.margin.electorate,elected,h.cv.elec,"Sao Paulo", bw.discont.samp = c(.006, .004, .002)))
occ.agricultor.check.elec  <-with(applied.data,discont.point.est(occ_agricultor,vote.margin.electorate,elected,h.cv.elec,"Occupation: Agriculture", bw.discont.samp = c(.006, .004, .002)))
occ.business.check.elec  <-with(applied.data,discont.point.est(occ_business,vote.margin.electorate,elected,h.cv.elec,"Occupation: Business", bw.discont.samp = c(.006, .004, .002)))
pmdb.check.elec  <-with(applied.data,discont.point.est(party_pmdb,vote.margin.electorate,elected,h.cv.elec,"Party: PMDB", bw.discont.samp = c(.006, .004, .002)))
psdb.check.elec  <-with(applied.data,discont.point.est(party_psdb,vote.margin.electorate,elected,h.cv.elec,"Party: PSDB", bw.discont.samp = c(.006, .004, .002)))
pfl.check.elec  <-with(applied.data,discont.point.est(party_pfl,vote.margin.electorate,elected,h.cv.elec,"Party: PFL", bw.discont.samp = c(.006, .004, .002)))
pt.pres.98.check.elec <-with(applied.data,discont.point.est(pt_pres_1998,vote.margin.electorate,elected,h.cv.elec,"Lula Vote Share, 1998", bw.discont.samp = c(.006, .004, .002)))
gini.check.elec  <-with(applied.data,discont.point.est(gini_2000,vote.margin.electorate,elected,h.cv.elec,"Gini (2000)", bw.discont.samp = c(.006, .004, .002)))
latitude.check.elec  <-with(applied.data,discont.point.est(latitude,vote.margin.electorate,elected,h.cv.elec,"Latitude", bw.discont.samp = c(.006, .004, .002)))
longitude.check.elec  <-with(applied.data,discont.point.est(longitude,vote.margin.electorate,elected,h.cv.elec,"Longitude", bw.discont.samp = c(.006, .004, .002)))
hdi.check.elec  <-with(applied.data,discont.point.est(hdi_2000,vote.margin.electorate,elected,h.cv.elec,"HDI (2000)", bw.discont.samp = c(.006, .004, .002)))
income.check.elec  <-with(applied.data,discont.point.est(income_2000,vote.margin.electorate,elected,h.cv.elec,"GDP per Cap (2000)", bw.discont.samp = c(.006, .004, .002)))
pt.pref.check.elec  <-with(applied.data,discont.point.est(pt_pref_2000,vote.margin.electorate,elected,h.cv.elec,"PT Mayor (2000)", bw.discont.samp = c(.006, .004, .002)))
time.since.app.check.elec <- with(applied.data,discont.point.est(time.since.app,vote.margin.electorate,elected,h.cv.elec,"Time Since Application", bw.discont.samp = c(.006, .004, .002)))


#Put together into one data.set and produce a latex table
bal.check.results.elec <- data.frame(rbind(apply.check.elec, prior.approval.check.elec,yob.check.elec,coalition.check.elec,edu.primary.check.elec,uf.sp.check.elec,uf.mg.check.elec,occ.agricultor.check.elec,occ.business.check.elec,pmdb.check.elec,psdb.check.elec,pfl.check.elec,pt.pres.98.check.elec,gini.check.elec,latitude.check.elec,longitude.check.elec,hdi.check.elec,pt.pref.check.elec,electorate.check.elec,time.since.app.check.elec))

n.var <- length(unique(bal.check.results.elec$variable))
bal.check.table.elec <- (matrix(nrow=n.var*2,ncol=6))
bal.check.table.elec[seq(1,n.var*2,by=2),1] <- unique(as.character(bal.check.results.elec$variable))
bal.check.table.elec[seq(2,n.var*2,by=2),1] <- ""
bal.check.table.elec[seq(1,n.var*2,by=2),2] <- bal.check.results.elec$coef[bal.check.results.elec$specification=="Full Sample"]
bal.check.table.elec[seq(1,n.var*2,by=2),3] <- bal.check.results.elec$coef[bal.check.results.elec$specification=="Local Linear Regression"]
bal.check.table.elec[seq(1,n.var*2,by=2),4] <- bal.check.results.elec$coef[bal.check.results.elec$specification=="Discontinuity Sample (vote margin < 0.006)"]
bal.check.table.elec[seq(1,n.var*2,by=2),5] <- bal.check.results.elec$coef[bal.check.results.elec$specification=="Discontinuity Sample (vote margin < 0.004)"]
bal.check.table.elec[seq(1,n.var*2,by=2),6] <- bal.check.results.elec$coef[bal.check.results.elec$specification=="Discontinuity Sample (vote margin < 0.002)"]
bal.check.table.elec[seq(2,n.var*2,by=2),2] <- bal.check.results.elec$se[bal.check.results.elec$specification=="Full Sample"]
bal.check.table.elec[seq(2,n.var*2,by=2),3] <- bal.check.results.elec$se[bal.check.results.elec$specification=="Local Linear Regression"]
bal.check.table.elec[seq(2,n.var*2,by=2),4] <- bal.check.results.elec$se[bal.check.results.elec$specification=="Discontinuity Sample (vote margin < 0.006)"]
bal.check.table.elec[seq(2,n.var*2,by=2),5] <- bal.check.results.elec$se[bal.check.results.elec$specification=="Discontinuity Sample (vote margin < 0.004)"]
bal.check.table.elec[seq(2,n.var*2,by=2),6] <- bal.check.results.elec$se[bal.check.results.elec$specification=="Discontinuity Sample (vote margin < 0.002)"]
colnames(bal.check.table.elec) <- c("Variable","Full Sample", "Loc. Lin. Reg", "Discont Sample (x<0.006)","Discont Sample (x<0.004)","Discont Sample (x<0.002)")
bal.check.table.elec[,2:6] <- apply(bal.check.table.elec[,2:6],2,function(x) round(as.numeric(x),digits=2))
bal.check.table.elec[seq(2,n.var*2,by=2),2:6] <- paste("(",bal.check.table.elec[seq(2,n.var*2,by=2),2:6],")",sep="")
bal.check.table.elec <- rbind(bal.check.table.elec[1:4,], c("n",unique(bal.check.results.elec$n[bal.check.results.elec$specification=="Full Sample"][1:2]),unique(bal.check.results.elec$n[bal.check.results.elec$specification=="Local Linear Regression"][1:2]),unique(bal.check.results.elec$n[bal.check.results.elec$specification=="Discontinuity Sample (vote margin < 0.006)"][1:2]),unique(bal.check.results.elec$n[bal.check.results.elec$specification=="Discontinuity Sample (vote margin < 0.004)"][1:2]),unique(bal.check.results.elec$n[bal.check.results.elec$specification=="Discontinuity Sample (vote margin < 0.002)"][1:2])), bal.check.table.elec[5:nrow(bal.check.table.elec),])
bal.check.table.elec[1,2:6] <- as.character(signif(apply.check.elec[c(5,1:4),3],2))
bal.check.table.elec[2,2:6] <- paste("(",as.character(signif(apply.check.elec[c(5, 1:4),4]),2),")",sep="")
bal.check.table.elec[3,2:6] <- as.character(signif(prior.approval.check.elec[c(5, 1:4),3],2))
bal.check.table.elec[4,2:6] <- paste("(",as.character(signif(prior.approval.check.elec[c(5, 1:4),4]),2),")",sep="")
bal.check.table.elec <- rbind(bal.check.table.elec, c("n",unique(bal.check.results.elec$n[bal.check.results.elec$specification=="Full Sample"][3]),unique(bal.check.results.elec$n[bal.check.results.elec$specification=="Local Linear Regression"][3]),unique(bal.check.results.elec$n[bal.check.results.elec$specification=="Discontinuity Sample (vote margin < 0.006)"][3]),unique(bal.check.results.elec$n[bal.check.results.elec$specification=="Discontinuity Sample (vote margin < 0.004)"][3]),unique(bal.check.results.elec$n[bal.check.results.elec$specification=="Discontinuity Sample (vote margin < 0.002)"][3])))

bal.check.table.elec <- xtable(bal.check.table.elec,display=c("s","s",rep("fg",5)), caption = "Balance Statistics")
digits(bal.check.table.elec) <-2
align(bal.check.table.elec) <- "rr|rrrrr"
print(bal.check.table.elec,latex.environment = "center", include.rownames = FALSE,hline.after=c(0,nrow(bal.check.table.elec)-1),floating.environment = "sidewaystable")


###
#Table 4
###
loc.lin.app.elec <- (lm(post.approved~elected*vote.margin.electorate,data=applied.data,subset=abs(vote.margin.electorate)<=h.cv.elec))
loc.lin.rej.elec <- (lm(1-pending_yes~elected*vote.margin.electorate,data=applied.data,subset=abs(vote.margin.electorate)<=h.cv.elec))
discont1.app.elec <- lm(post.approved~elected, data=applied.data,subset = abs(vote.margin.electorate)<=.006)
discont1.rej.elec <- lm(1-pending_yes~elected, data=applied.data,subset = abs(vote.margin.electorate)<=.006)
discont2.app.elec <- lm(post.approved~elected, data=applied.data,subset = abs(vote.margin.electorate)<=.004)
discont2.rej.elec <- lm(1-pending_yes~elected, data=applied.data,subset = abs(vote.margin.electorate)<=0.004)
discont3.app.elec <- lm(post.approved~elected, data=applied.data,subset = abs(vote.margin.electorate)<=.002)
discont3.rej.elec <- lm(1-pending_yes~elected, data=applied.data,subset = abs(vote.margin.electorate)<=.002)

app.coefs.elec <- c(coef(loc.lin.app.elec)[2],coef(discont1.app.elec)[2],coef(discont2.app.elec)[2],coef(discont3.app.elec)[2])
app.se.elec <- sqrt(c(vcov(loc.lin.app.elec)[2,2],vcov(discont1.app.elec)[2,2],vcov(discont2.app.elec)[2,2],vcov(discont3.app.elec)[2,2]))
app.n.elec <- c(nrow(loc.lin.app.elec$model),nrow(discont1.app.elec$model),nrow(discont2.app.elec$model), nrow(discont3.app.elec$model))
app.table.elec <- (matrix(rbind(app.coefs.elec, app.se.elec, app.n.elec),nrow=3))

rej.coefs.elec <- c(coef(loc.lin.rej.elec)[2],coef(discont1.rej.elec)[2],coef(discont2.rej.elec)[2],coef(discont1.rej.elec)[2])
rej.se.elec <- sqrt(c(vcov(loc.lin.rej.elec)[2,2],vcov(discont1.rej.elec)[2,2],vcov(discont2.rej.elec)[2,2],vcov(discont3.rej.elec)[2,2]))
rej.n.elec <- c(nrow(loc.lin.rej.elec$model),nrow(discont1.rej.elec$model),nrow(discont2.rej.elec$model), nrow(discont3.rej.elec$model))
rej.table.elec <- (matrix(rbind(rej.coefs.elec,rej.se.elec,rej.n.elec),nrow=3))

xtable(rbind(app.table.elec, rej.table.elec))

##
#Table 5
##
#Calculate inflation factor and inflated vote margin
full.data$infl.fact<-(full.data$votes/(full.data$votes-full.data$raw.vote.margin))
full.data$infl.fact[full.data$elected==0] <- full.data$infl.fact[full.data$elected==0]^-1 
full.data$infl.vote.margin <- full.data$raw.vote.margin * full.data$infl.fact
applied.data$infl.fact<-(applied.data$votes/(applied.data$votes-applied.data$raw.vote.margin))
applied.data$infl.fact[applied.data$elected==0]<-applied.data$infl.fact[applied.data$elected==0]^-1 
applied.data$infl.vote.margin<-applied.data$raw.vote.margin*applied.data$infl.fact

##Results###
approved.results <- with(applied.data,discont.point.est(y = post.approved, force.var = infl.vote.margin,treat = elected,bw.loc.lin= h.cv, y.name ="License Approval"))
rejected.results <- with(applied.data,discont.point.est(y = 1-pending_yes, force.var = infl.vote.margin,treat = elected,bw.loc.lin= h.cv, y.name ="License Rejection"))

##Covariates##
#Decision to apply
apply.check.infl <- with(full.data, discont.point.est(y = post.applied, force.var = infl.vote.margin, treat = elected, bw.loc.lin = h.cv, y.name = "Apply"))
#Prior approval
prior.approval.check.infl <-  with(full.data, discont.point.est(y = pre.approved, force.var = infl.vote.margin, treat = elected, bw.loc.lin = h.cv, y.name = "Prior Approval"))
#other covariates
yob.check.infl <- with(applied.data,discont.point.est(yob,infl.vote.margin,elected,h.cv,"Year of Birth"))
electorate.check.infl <- with(applied.data,discont.point.est(log(total.votes),infl.vote.margin,elected,h.cv,"Log Electorate"))
coalition.check.infl <- with(applied.data,discont.point.est(log(coal.votes),infl.vote.margin,elected,h.cv,"Log Coalition Votes"))
edu.primary.check.infl <- with(applied.data,discont.point.est(edu_primary,infl.vote.margin,elected,h.cv,"Primary Education"))
pre.approved.check.infl <- with(applied.data,discont.point.est(pre.approved,infl.vote.margin,elected,h.cv,"Approved Before Election"))
uf.mg.check.infl <- with(applied.data,discont.point.est(uf_mg,infl.vote.margin,elected,h.cv,"Minas Gerais"))
uf.sp.check.infl <- with(applied.data,discont.point.est(uf_sp,infl.vote.margin,elected,h.cv,"Sao Paulo"))
occ.agricultor.check.infl  <-with(applied.data,discont.point.est(occ_agricultor,infl.vote.margin,elected,h.cv,"Occupation: Agriculture"))
occ.business.check.infl  <-with(applied.data,discont.point.est(occ_business,infl.vote.margin,elected,h.cv,"Occupation: Business"))
pmdb.check.infl  <-with(applied.data,discont.point.est(party_pmdb,infl.vote.margin,elected,h.cv,"Party: PMDB"))
psdb.check.infl  <-with(applied.data,discont.point.est(party_psdb,infl.vote.margin,elected,h.cv,"Party: PSDB"))
pfl.check.infl  <-with(applied.data,discont.point.est(party_pfl,infl.vote.margin,elected,h.cv,"Party: PFL"))
pt.pres.98.check.infl <-with(applied.data,discont.point.est(pt_pres_1998,infl.vote.margin,elected,h.cv,"Lula Vote Share, 1998"))
gini.check.infl  <-with(applied.data,discont.point.est(gini_2000,infl.vote.margin,elected,h.cv,"Gini (2000)"))
latitude.check.infl  <-with(applied.data,discont.point.est(latitude,infl.vote.margin,elected,h.cv,"Latitude"))
longitude.check.infl  <-with(applied.data,discont.point.est(longitude,infl.vote.margin,elected,h.cv,"Longitude"))
hdi.check.infl  <-with(applied.data,discont.point.est(hdi_2000,infl.vote.margin,elected,h.cv,"HDI (2000)"))
income.check.infl  <-with(applied.data,discont.point.est(income_2000,infl.vote.margin,elected,h.cv,"GDP per Cap (2000)"))
pt.pref.check.infl  <-with(applied.data,discont.point.est(pt_pref_2000,infl.vote.margin,elected,h.cv,"PT Mayor (2000)"))
time.since.app.check.infl <- with(applied.data,discont.point.est(time.since.app,infl.vote.margin,elected,h.cv,"Time Since Application"))


#Put together into one data.set and produce a latex table
bal.check.results.infl <- data.frame(rbind(apply.check.infl, prior.approval.check.infl,yob.check.infl,coalition.check.infl,edu.primary.check.infl,uf.sp.check.infl,uf.mg.check.infl,occ.agricultor.check.infl,occ.business.check.infl,pmdb.check.infl,psdb.check.infl,pfl.check.infl,pt.pres.98.check.infl,gini.check.infl,latitude.check.infl,longitude.check.infl,hdi.check.infl,pt.pref.check.infl,electorate.check.infl,time.since.app.check.infl))

n.var <- length(unique(bal.check.results.infl$variable))
bal.check.table.infl <- (matrix(nrow=n.var*2,ncol=6))
bal.check.table.infl[seq(1,n.var*2,by=2),1] <- unique(as.character(bal.check.results.infl$variable))
bal.check.table.infl[seq(2,n.var*2,by=2),1] <- ""
bal.check.table.infl[seq(1,n.var*2,by=2),2] <- bal.check.results.infl$coef[bal.check.results.infl$specification=="Full Sample"]
bal.check.table.infl[seq(1,n.var*2,by=2),3] <- bal.check.results.infl$coef[bal.check.results.infl$specification=="Local Linear Regression"]
bal.check.table.infl[seq(1,n.var*2,by=2),4] <- bal.check.results.infl$coef[bal.check.results.infl$specification=="Discontinuity Sample (vote margin < 40)"]
bal.check.table.infl[seq(1,n.var*2,by=2),5] <- bal.check.results.infl$coef[bal.check.results.infl$specification=="Discontinuity Sample (vote margin < 20)"]
bal.check.table.infl[seq(1,n.var*2,by=2),6] <- bal.check.results.infl$coef[bal.check.results.infl$specification=="Discontinuity Sample (vote margin < 10)"]
bal.check.table.infl[seq(2,n.var*2,by=2),2] <- bal.check.results.infl$se[bal.check.results.infl$specification=="Full Sample"]
bal.check.table.infl[seq(2,n.var*2,by=2),3] <- bal.check.results.infl$se[bal.check.results.infl$specification=="Local Linear Regression"]
bal.check.table.infl[seq(2,n.var*2,by=2),4] <- bal.check.results.infl$se[bal.check.results.infl$specification=="Discontinuity Sample (vote margin < 10)"]
bal.check.table.infl[seq(2,n.var*2,by=2),5] <- bal.check.results.infl$se[bal.check.results.infl$specification=="Discontinuity Sample (vote margin < 20)"]
bal.check.table.infl[seq(2,n.var*2,by=2),6] <- bal.check.results.infl$se[bal.check.results.infl$specification=="Discontinuity Sample (vote margin < 40)"]
colnames(bal.check.table.infl) <- c("Variable","Full Sample", "Loc. Lin. Reg", "Discont Sample (x<40)","Discont Sample (x<20)","Discont Sample (x<10)")
bal.check.table.infl[,2:6] <- apply(bal.check.table.infl[,2:6],2,function(x) round(as.numeric(x),digits=2))
bal.check.table.infl[seq(2,n.var*2,by=2),2:6] <- paste("(",bal.check.table.infl[seq(2,n.var*2,by=2),2:6],")",sep="")

bal.check.table.infl <- rbind(bal.check.table.infl[1:4,], c("n",unique(bal.check.results.infl$n[bal.check.results.infl$specification=="Full Sample"][1:2]),unique(bal.check.results.infl$n[bal.check.results.infl$specification=="Local Linear Regression"][1:2]),unique(bal.check.results.infl$n[bal.check.results.infl$specification=="Discontinuity Sample (vote margin < 40)"][1:2]),unique(bal.check.results.infl$n[bal.check.results.infl$specification=="Discontinuity Sample (vote margin < 20)"][1:2]),unique(bal.check.results.infl$n[bal.check.results.infl$specification=="Discontinuity Sample (vote margin < 10)"][1:2])), bal.check.table.infl[5:nrow(bal.check.table.infl),])
bal.check.table.infl[1,2:6] <- as.character(signif(apply.check.infl[c(5,1:4),3],2))
bal.check.table.infl[2,2:6] <- paste("(",as.character(signif(apply.check.infl[c(5,1:4),4]),2),")",sep="")
bal.check.table.infl[3,2:6] <- as.character(signif(prior.approval.check.infl[c(5,1:4),3],2))
bal.check.table.infl[4,2:6] <- paste("(",as.character(signif(prior.approval.check.infl[c(5,1:4),4]),2),")",sep="")
bal.check.table.infl <- rbind(bal.check.table.infl, c("n",unique(bal.check.results.infl$n[bal.check.results.infl$specification=="Full Sample"][3]),unique(bal.check.results.infl$n[bal.check.results.infl$specification=="Local Linear Regression"][3]),unique(bal.check.results.infl$n[bal.check.results.infl$specification=="Discontinuity Sample (vote margin < 40)"][3]),unique(bal.check.results.infl$n[bal.check.results.infl$specification=="Discontinuity Sample (vote margin < 20)"][3]),unique(bal.check.results.infl$n[bal.check.results.infl$specification=="Discontinuity Sample (vote margin < 10)"][3])))



bal.check.table.infl <- xtable(bal.check.table.infl,display=c("s","s",rep("fg",5)), caption = "Balance Statistics")
digits(bal.check.table.infl) <-2
align(bal.check.table.infl) <- "rr|rrrrr"
print(bal.check.table.infl,latex.environment = "center", include.rownames = FALSE,hline.after=c(0,nrow(bal.check.table.infl)-1),floating.environment = "sidewaystable")


##
#Table 6
##
loc.lin.app.infl <- (lm(post.approved~elected*infl.vote.margin,data=applied.data,subset=abs(infl.vote.margin)<=h.cv))
loc.lin.rej.infl<- (lm(1-pending_yes~elected*infl.vote.margin,data=applied.data,subset=abs(infl.vote.margin)<=h.cv))
discont40.app.infl <- lm(post.approved~elected, data=applied.data,subset = abs(infl.vote.margin)<=40)
discont40.rej.infl <- lm(1-pending_yes~elected, data=applied.data,subset = abs(infl.vote.margin)<=40)
discont20.app.infl <- lm(post.approved~elected, data=applied.data,subset = abs(infl.vote.margin)<=20)
discont20.rej.infl <- lm(1-pending_yes~elected, data=applied.data,subset = abs(infl.vote.margin)<=20)
discont10.app.infl <- lm(post.approved~elected, data=applied.data,subset = abs(infl.vote.margin)<=10)
discont10.rej.infl <- lm(1-pending_yes~elected, data=applied.data,subset = abs(infl.vote.margin)<=10)

app.coefs <- c(coef(loc.lin.app.infl)[2],coef(discont40.app.infl)[2],coef(discont20.app.infl)[2],coef(discont10.app.infl)[2])
app.se <- sqrt(c(vcov(loc.lin.app.infl)[2,2],vcov(discont40.app.infl)[2,2],vcov(discont20.app.infl)[2,2],vcov(discont10.app.infl)[2,2]))
app.n <- c(nrow(loc.lin.app.infl$model),nrow(discont40.app.infl$model),nrow(discont20.app.infl$model), nrow(discont10.app.infl$model))
app.table <- (matrix(rbind(app.coefs,app.se,app.n),nrow=3))
rej.coefs <- c(coef(loc.lin.rej.infl)[2],coef(discont40.rej.infl)[2],coef(discont20.rej.infl)[2],coef(discont10.rej.infl)[2])
rej.se <- sqrt(c(vcov(loc.lin.rej.infl)[2,2],vcov(discont40.rej.infl)[2,2],vcov(discont20.rej.infl)[2,2],vcov(discont10.rej.infl)[2,2]))
rej.n <- c(nrow(loc.lin.rej.infl$model),nrow(discont40.rej.infl$model),nrow(discont20.rej.infl$model), nrow(discont10.rej.infl$model))
rej.table <- (matrix(rbind(rej.coefs,rej.se,rej.n),nrow=3))
xtable(rbind(app.table,rej.table))

###
##Table 7
##

loc.lin.small <- (lm(post.approved~elected*vote.margin.electorate,data=applied.data[applied.data$pop_2000<14149,],subset=abs(vote.margin.electorate)<=h.cv.elec))
loc.lin.big <- (lm(post.approved~elected*vote.margin.electorate,data=applied.data[applied.data$pop_2000>14149,],subset=abs(vote.margin.electorate)<=h.cv.elec))
discont1.small <- lm(post.approved~elected, data=applied.data[applied.data$pop_2000<14149,],subset = abs(vote.margin.electorate)<=0.006)
discont1.big <- lm(post.approved~elected, data=applied.data[applied.data$pop_2000>14149,],subset = abs(vote.margin.electorate)<=0.006)
discont2.small <- lm(post.approved~elected, data=applied.data[applied.data$pop_2000<14149,],subset = abs(vote.margin.electorate)<=0.004)
discont2.big <- lm(post.approved~elected, data=applied.data[applied.data$pop_2000>14149,],subset = abs(vote.margin.electorate)<=0.004)
discont3.small <- lm(post.approved~elected, data=applied.data[applied.data$pop_2000<14149,],subset = abs(vote.margin.electorate)<=0.002)
discont3.big <- lm(post.approved~elected, data=applied.data[applied.data$pop_2000>14149,],subset = abs(vote.margin.electorate)<=0.002)

small.coefs <- c(coef(loc.lin.small)[2],coef(discont1.small)[2],coef(discont2.small)[2],coef(discont3.small)[2])
small.se <- sqrt(c(vcov(loc.lin.small)[2,2],vcov(discont1.small)[2,2],vcov(discont2.small)[2,2],vcov(discont3.small)[2,2]))
small.n <- c(nrow(loc.lin.small$model),nrow(discont1.small$model),nrow(discont2.small$model), nrow(discont3.small$model))
small.table.elec <- (matrix(rbind(small.coefs,small.se,small.n),nrow=3))

big.coefs <- c(coef(loc.lin.big)[2],coef(discont1.big)[2],coef(discont2.big)[2],coef(discont3.big)[2])
big.se <- sqrt(c(vcov(loc.lin.big)[2,2],vcov(discont1.big)[2,2],vcov(discont2.big)[2,2],vcov(discont3.big)[2,2]))
big.n <- c(nrow(loc.lin.big$model),nrow(discont1.big$model),nrow(discont2.big$model), nrow(discont3.big$model))
big.table.elec <- (matrix(rbind(big.coefs,big.se,big.n),nrow=3))

loc.lin.small.rej <- (lm(1-pending_yes~elected*vote.margin.electorate,data=applied.data[applied.data$pop_2000<14149,],subset=abs(vote.margin.electorate)<=h.cv.elec))
loc.lin.big.rej <- (lm(1-pending_yes~elected*vote.margin.electorate,data=applied.data[applied.data$pop_2000>14149,],subset=abs(vote.margin.electorate)<=h.cv.elec))
discont1.small.rej <- lm(1-pending_yes~elected, data=applied.data[applied.data$pop_2000<14149,],subset = abs(vote.margin.electorate)<=0.006)
discont1.big.rej <- lm(1-pending_yes~elected, data=applied.data[applied.data$pop_2000>14149,],subset = abs(vote.margin.electorate)<=0.006)
discont2.small.rej <- lm(1-pending_yes~elected, data=applied.data[applied.data$pop_2000<14149,],subset = abs(vote.margin.electorate)<=0.004)
discont2.big.rej <- lm(1-pending_yes~elected, data=applied.data[applied.data$pop_2000>14149,],subset = abs(vote.margin.electorate)<=0.004)
discont3.small.rej <- lm(1-pending_yes~elected, data=applied.data[applied.data$pop_2000<14149,],subset = abs(vote.margin.electorate)<=0.002)
discont3.big.rej <- lm(1-pending_yes~elected, data=applied.data[applied.data$pop_2000>14149,],subset = abs(vote.margin.electorate)<=0.002)

small.coefs <- c(coef(loc.lin.small.rej)[2],coef(discont1.small.rej)[2],coef(discont2.small.rej)[2],coef(discont3.small.rej)[2])
small.se <- sqrt(c(vcov(loc.lin.small.rej)[2,2],vcov(discont1.small.rej)[2,2],vcov(discont2.small.rej)[2,2],vcov(discont3.small.rej)[2,2]))
small.n <- c(nrow(loc.lin.small.rej$model),nrow(discont1.small.rej$model),nrow(discont2.small.rej$model), nrow(discont3.small.rej$model))
small.table.rej.elec <- (matrix(rbind(small.coefs,small.se,small.n),nrow=3))

big.coefs <- c(coef(loc.lin.big.rej)[2],coef(discont1.big.rej)[2],coef(discont2.big.rej)[2],coef(discont3.big.rej)[2])
big.se <- sqrt(c(vcov(loc.lin.big.rej)[2,2],vcov(discont1.big.rej)[2,2],vcov(discont2.big.rej)[2,2],vcov(discont3.big.rej)[2,2]))
big.n <- c(nrow(loc.lin.big.rej$model),nrow(discont1.big.rej$model),nrow(discont2.big.rej$model), nrow(discont3.big.rej$model))
big.table.rej.elec <- (matrix(rbind(big.coefs,big.se,big.n),nrow=3))
xtable(rbind(small.table.elec,big.table.elec, small.table.rej.elec, big.table.rej.elec))


#
#Table 8
#

##Results###
loc.lin.small <- (lm(post.approved~elected*raw.vote.margin,data=applied.data[applied.data$pop_2000<14149,],subset=abs(raw.vote.margin)<=h.cv))
loc.lin.big <- (lm(post.approved~elected*raw.vote.margin,data=applied.data[applied.data$pop_2000>14149,],subset=abs(raw.vote.margin)<=h.cv))
discont40.small <- lm(post.approved~elected, data=applied.data[applied.data$pop_2000<14149,],subset = abs(raw.vote.margin)<=40)
discont40.big <- lm(post.approved~elected, data=applied.data[applied.data$pop_2000>14149,],subset = abs(raw.vote.margin)<=40)
discont20.small <- lm(post.approved~elected, data=applied.data[applied.data$pop_2000<14149,],subset = abs(raw.vote.margin)<=20)
discont20.big <- lm(post.approved~elected, data=applied.data[applied.data$pop_2000>14149,],subset = abs(raw.vote.margin)<=20)
discont10.small <- lm(post.approved~elected, data=applied.data[applied.data$pop_2000<14149,],subset = abs(raw.vote.margin)<=10)
discont10.big <- lm(post.approved~elected, data=applied.data[applied.data$pop_2000>14149,],subset = abs(raw.vote.margin)<=10)

small.coefs <- c(coef(loc.lin.small)[2],coef(discont40.small)[2],coef(discont20.small)[2],coef(discont10.small)[2])
small.se <- sqrt(c(vcov(loc.lin.small)[2,2],vcov(discont40.small)[2,2],vcov(discont20.small)[2,2],vcov(discont10.small)[2,2]))
small.n <- c(nrow(loc.lin.small$model),nrow(discont40.small$model),nrow(discont20.small$model), nrow(discont10.small$model))
small.table <- (matrix(rbind(small.coefs,small.se,small.n),nrow=3))

big.coefs <- c(coef(loc.lin.big)[2],coef(discont40.big)[2],coef(discont20.big)[2],coef(discont10.big)[2])
big.se <- sqrt(c(vcov(loc.lin.big)[2,2],vcov(discont40.big)[2,2],vcov(discont20.big)[2,2],vcov(discont10.big)[2,2]))
big.n <- c(nrow(loc.lin.big$model),nrow(discont40.big$model),nrow(discont20.big$model), nrow(discont10.big$model))
big.table <- (matrix(rbind(big.coefs,big.se,big.n),nrow=3))
xtable(rbind(small.table,big.table))

#Application rejection
##Results###
loc.lin.small.rej <- (lm(1-pending_yes~elected*raw.vote.margin,data=applied.data[applied.data$pop_2000<14149,],subset=abs(raw.vote.margin)<=h.cv))
loc.lin.big.rej <- (lm(1-pending_yes~elected*raw.vote.margin,data=applied.data[applied.data$pop_2000>14149,],subset=abs(raw.vote.margin)<=h.cv))
discont40.small.rej <- lm(1-pending_yes~elected, data=applied.data[applied.data$pop_2000<14149,],subset = abs(raw.vote.margin)<=40)
discont40.big.rej <- lm(1-pending_yes~elected, data=applied.data[applied.data$pop_2000>14149,],subset = abs(raw.vote.margin)<=40)
discont20.small.rej <- lm(1-pending_yes~elected, data=applied.data[applied.data$pop_2000<14149,],subset = abs(raw.vote.margin)<=20)
discont20.big.rej <- lm(1-pending_yes~elected, data=applied.data[applied.data$pop_2000>14149,],subset = abs(raw.vote.margin)<=20)
discont10.small.rej <- lm(1-pending_yes~elected, data=applied.data[applied.data$pop_2000<14149,],subset = abs(raw.vote.margin)<=10)
discont10.big.rej <- lm(1-pending_yes~elected, data=applied.data[applied.data$pop_2000>14149,],subset = abs(raw.vote.margin)<=10)

small.coefs <- c(coef(loc.lin.small.rej)[2],coef(discont40.small.rej)[2],coef(discont20.small.rej)[2],coef(discont10.small.rej)[2])
small.se <- sqrt(c(vcov(loc.lin.small.rej)[2,2],vcov(discont40.small.rej)[2,2],vcov(discont20.small.rej)[2,2],vcov(discont10.small.rej)[2,2]))
small.n <- c(nrow(loc.lin.small.rej$model),nrow(discont40.small.rej$model),nrow(discont20.small.rej$model), nrow(discont10.small.rej$model))
small.table.rej <- (matrix(rbind(small.coefs,small.se,small.n),nrow=3))

big.coefs <- c(coef(loc.lin.big.rej)[2],coef(discont40.big.rej)[2],coef(discont20.big.rej)[2],coef(discont10.big.rej)[2])
big.se <- sqrt(c(vcov(loc.lin.big.rej)[2,2],vcov(discont40.big.rej)[2,2],vcov(discont20.big.rej)[2,2],vcov(discont10.big.rej)[2,2]))
big.n <- c(nrow(loc.lin.big.rej$model),nrow(discont40.big.rej$model),nrow(discont20.big.rej$model), nrow(discont10.big.rej$model))
big.table.rej <- (matrix(rbind(big.coefs,big.se,big.n),nrow=3))

xtable(rbind(small.table, big.table, small.table.rej, big.table.rej))



###
#Table 9
###
match.pctVV.adj <- Match(Y= match.data$pctVV, Tr = match.data$treat, X = genmatch.covar, Weight.matrix=genmatch.output, Z=data.frame( match.data$log.total.assets), BiasAdjust=TRUE,  caliper = c(rep(10,43),.5))
match.elected.adj <-  Match(Y= match.data$elected, Tr = match.data$treat, X = genmatch.covar, Weight.matrix=genmatch.output, Z=data.frame( match.data$log.total.assets), BiasAdjust=TRUE,  caliper = c(rep(10,43),.5))

results.tab.adj <- data.frame("Pct Valid Votes" = c((match.pctVV.adj$est), (match.pctVV.adj$se), (length(match.pctVV.adj$mdata$Y))),
"Elected" = c((match.elected.adj$est), (match.elected.adj$se), (length(match.elected.adj$mdata$Y))))
rownames(results.tab.adj) <- c( "ATT Estimate", "SE", "n")

caliper <- c(100, .5, 100, 100, 1, 1 , 100, 100, 1.5, 100, 100, 100, 100, 100, 100, 100,100)
exact <- (c(0,0,1,0,0,1,0,0,0,0,0,0,0,1,1,0,0))

match.pctVV.alt1 <- Match(Y=match.data.alt1$pctVV, Tr=match.data.alt1$treat, X= genmatch.covar.alt1 ,M=1, estimand="ATT", Weight.matrix=genmatch.output.alt1, caliper = caliper, exact = exact)
match.elected.alt1 <-  Match(Y=match.data.alt1$elected, Tr=match.data.alt1$treat, X= genmatch.covar.alt1 ,M=1, estimand="ATT", Weight.matrix=genmatch.output.alt1, caliper = caliper, exact = exact)

results.tab.alt1 <- data.frame("Pct Valid Votes" = c((match.pctVV.alt1$est), (match.pctVV.alt1$se), (length(match.pctVV.alt1$mdata$Y))), "Elected" = c((match.elected.alt1$est), (match.elected.alt1$se), (length(match.elected.alt1$mdata$Y))))
rownames(results.tab.alt1) <- c( "ATT Estimate", "SE", "n")



match.pctVV.alt2 <- Match(Y=match.data.alt2$pctVV, Tr=match.data.alt2$treat, X= genmatch.covar.alt2, M=1, estimand="ATT", Weight.matrix=genmatch.output.alt2)
match.elected.alt2 <- Match(Y=match.data.alt2$elected, Tr=match.data.alt2$treat, X= genmatch.covar.alt2, M=1, estimand="ATT", Weight.matrix=genmatch.output.alt2)

results.tab.alt2 <- data.frame("Pct Valid Votes" = c((match.pctVV.alt2$est), (match.pctVV.alt2$se), (length(match.pctVV.alt2$mdata$Y))), "Elected" = c((match.elected.alt2$est), (match.elected.alt2$se), (length(match.elected.alt2$mdata$Y))))
rownames(results.tab.alt2) <- c( "ATT Estimate", "SE", "n")


results.tab.alt <-  rbind(results.tab.adj, results.tab.alt1, results.tab.alt2)
results.tab.alt <- xtable( results.tab.alt, display=c("s",rep("fg",2)), caption = "Estimates")
digits(results.tab.alt) <-2
align(results.tab.alt) <- "r|rr"
print(results.tab.alt,latex.environment = "center", include.rownames = TRUE,hline.after=c(0))



##
#Table 10
##
bal.data <- with(match.data.alt1, data.frame(party, male, occ, edu, yob, lat, long, ran.prior, incumbent, log.valid.votes, prior.pctVV, elec.year, uf.ba, uf.sp, uf.mg, uf.rs , party.prior.pctVV, gini_2000, pt_pres_1998, income_2000, psdb_2000, pt_2000, hdi_2000, log.total.assets, time.since.app, log.num.apps))
bal.fmla <- as.formula(paste("treat~",paste(names(bal.data),collapse="+")))
match.alt1.bal <- MatchBalance( bal.fmla, match.out=match.pctVV.alt1, data=match.data.alt1, nboots=3000)

before.sd.diff <- sapply(match.alt1.bal$BeforeMatching, function(x) x[[1]])
after.sd.diff <- sapply(match.alt1.bal$AfterMatching, function(x) x[[1]])
before.t.p <- matrix(sapply(match.alt1.bal$BeforeMatching, function(x) x$tt[3]))
after.t.p <- matrix(sapply(match.alt1.bal$AfterMatching, function(x) x$tt[3]))
before.ks.p <- matrix(sapply(match.alt1.bal$BeforeMatching, function(x) x$ks$ks.boot.pvalue))
after.ks.p <- matrix(sapply(match.alt1.bal$AfterMatching, function(x) x$ks$ks.boot.pvalue))

bal.stats <- data.frame(var = names(data.frame(model.matrix(as.formula(paste("~",paste(names(bal.data),collapse="+"))), bal.data)[,-1])), before.sd.diff, after.sd.diff, before.t.p = unlist(before.t.p), after.t.p = unlist(after.t.p), before.ks.p, after.ks.p)
bal.stats <- bal.stats[c(4,7,18,21,27:31,33:36, 38:45, 51, 46, 58, 52:57,47:50, 60),]
rownames(bal.stats) <- c("Party: PFL", "Party: PMDB", "Party: PSDB", "Party: PT", "Male", "Occupation: Blue Collar", "Occupation: Education", "Occupation: Government", "Occupation: Media", "Occupation: None", "Occupation: Other", "Occupation: Politician", "Occupation: White Collar", "Education: Some Superior or More", "Year of Birth", "Latitude", "Longitude", "Ran Previously", "Incumbency", "Log Electorate", "Prior Vote Share", "Party Prior Vote Share", "Election Year", "Total Assets", "2000 Gini", "PT Pres Vote Share (1998)", "GDP Per Capita (2000)", "PSDB Mayor Vote Share (2000)", "PT Mayor Vote Share (2000)", "HDI (2000)", "State: Bahia", "State: Sao Paulo", "State: Minas Gerais", "State: Rio Grande do Sul", "Log Number of Applications")
bal.stats$before.ks.p[bal.stats$before.ks.p=="NULL"] <- NA
bal.stats$after.ks.p[bal.stats$after.ks.p=="NULL"] <- NA
bal.stats$before.ks.p <- unlist(bal.stats$before.ks.p)
bal.stats$after.ks.p <- unlist(bal.stats$after.ks.p)
bal.stats <- bal.stats[order(abs(bal.stats$before.sd.diff), decreasing = TRUE),]
bal.stats <- bal.stats[,-1]

bal.check.table <- xtable(bal.stats, display=c("s",rep("fg",6)), caption = "Balance Statistics")
digits(bal.check.table) <-2
align(bal.check.table) <- "r|rrrrrr"
print(bal.check.table,latex.environment = "center", include.rownames = TRUE,hline.after=c(0),floating.environment = "sidewaystable")


###
#Table 11
###

bal.data <- with(match.data.alt2, data.frame(party, male, occ, edu, yob, lat, long, ran.prior, incumbent, log.valid.votes, prior.pctVV, elec.year, uf.ba, uf.sp, uf.mg, uf.rs , party.prior.pctVV, gini_2000, pt_pres_1998, income_2000, psdb_2000, pt_2000, hdi_2000, log.total.assets, time.since.app, log.num.apps))
bal.fmla <- as.formula(paste("treat~",paste(names(bal.data),collapse="+")))
matchalt2.bal <- MatchBalance( bal.fmla, match.out=match.pctVV.alt2, data=match.data.alt2, nboots=3000)

before.sd.diff <- sapply(matchalt2.bal$BeforeMatching, function(x) x[[1]])
after.sd.diff <- sapply(matchalt2.bal$AfterMatching, function(x) x[[1]])
before.t.p <- matrix(sapply(matchalt2.bal$BeforeMatching, function(x) x$tt[3]))
after.t.p <- matrix(sapply(matchalt2.bal$AfterMatching, function(x) x$tt[3]))
before.ks.p <- matrix(sapply(matchalt2.bal$BeforeMatching, function(x) x$ks$ks.boot.pvalue))
after.ks.p <- matrix(sapply(matchalt2.bal$AfterMatching, function(x) x$ks$ks.boot.pvalue))

bal.stats <- data.frame(var = names(data.frame(model.matrix(as.formula(paste("~",paste(names(bal.data),collapse="+"))), bal.data)[,-1])), before.sd.diff, after.sd.diff, before.t.p = unlist(before.t.p), after.t.p = unlist(after.t.p), before.ks.p, after.ks.p)
bal.stats <- bal.stats[c(4,7,18,22,28,30:37, 39:46, 52, 47, 59, 53:58,48:51, 61),]
rownames(bal.stats) <- c("Party: PFL", "Party: PMDB", "Party: PSDB", "Party: PT", "Male", "Occupation: Blue Collar", "Occupation: Education", "Occupation: Government", "Occupation: Media", "Occupation: None", "Occupation: Other", "Occupation: Politician", "Occupation: White Collar", "Education: Some Superior or More", "Year of Birth", "Latitude", "Longitude", "Ran Previously", "Incumbency", "Log Electorate", "Prior Vote Share", "Party Prior Vote Share", "Election Year", "Total Assets", "2000 Gini", "PT Pres Vote Share (1998)", "GDP Per Capita (2000)", "PSDB Mayor Vote Share (2000)", "PT Mayor Vote Share (2000)", "HDI (2000)", "State: Bahia", "State: Sao Paulo", "State: Minas Gerais", "State: Rio Grande do Sul", "Log Number of Applications")
bal.stats$before.ks.p[bal.stats$before.ks.p=="NULL"] <- NA
bal.stats$after.ks.p[bal.stats$after.ks.p=="NULL"] <- NA
bal.stats$before.ks.p <- unlist(bal.stats$before.ks.p)
bal.stats$after.ks.p <- unlist(bal.stats$after.ks.p)
bal.stats <- bal.stats[order(abs(bal.stats$before.sd.diff), decreasing = TRUE),]
bal.stats <- bal.stats[,-1]

bal.check.table <- xtable(bal.stats, display=c("s",rep("fg",6)), caption = "Balance Statistics")
digits(bal.check.table) <-2
align(bal.check.table) <- "r|rrrrrr"
print(bal.check.table,latex.environment = "center", include.rownames = TRUE,hline.after=c(0),floating.environment = "sidewaystable")

