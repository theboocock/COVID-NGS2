# This script plots the epidemiological parameters of a BDSKY analysis with their HPD intervals. 
# Written by Denise Kühnert (denise.kuehnert@gmail.com)
#
# start R from Terminal and run script with commandline: source(bdsky_plot.R) 
#
# Input: 
#	- a text file listing the log files that should be analyzed (see example file loglist.txt)
#	- burnin (the percentage of samples that should be ignored from the log file)
#	- date of the most recent sample (for plotting from past to present)
#
# Assumptions:
# - The given log file contains parameters named "reproductiveNumber" or "reproductiveNumberx" with index x,
#	"becomeUninfectiousRate" / "becomeUninfectiousRatex" and 
#	"samplingProportion" / "samplingProportionx"
# - The number of reproductiveNumber parameters determines the intervalnumber
# - The parameters becomeUninfectiousRate and samplingProportion are assumed to be either constant or to have the same dimension as reproductiveNumber.

a= 0.8e-3 
b =5e-4
scale =  b^2/a
k = a / scale 

#(a^2/b^2)
#sqrt(b^2/a
library(s20x)
library(boa)
library(Hmisc)
library(miscTools)
# input : a matrix M and a column ascii name 
# output : the numeric index of the column 
colnameindex = function(M , colname0) 
{ 
  colsnames = names(M[1,]); 
  theindex = which(colsnames==colname0); 
  return(theindex); 
}

# /* get user input */
#cat("Please enter the name of the file that stores the path to all log files to be analyzed: ")
#loglist <- data.frame(in_f=c("data/beast_log/b143_no_estimate_infectious_period.log","data/beast_phylo/la.log","data/strains_sub.log","data/beast_phylo/test_all_la.log"))

#cat("Please enter the percentage you would like to use as burnin (e.g. 10 for 10%): ")
burninpercent <- 10

#cat("Please enter the sampling date of the most recent sample: ")
recent <- "2020-06-28"

#cat("Please enter the grid size (can be bigger than the intervalnumber to get a smoother plot): ")
gridSize <- 20

# /* read log list */ 
#loglist = read.table(x, as.is=TRUE, header=FALSE)

#closeAllConnections()

process_beast_phylo = function(beast_phylo, recent,plot=T){
  
  # recent <- "2020-06-28"
  
  #bdsky_plot = function(log)
  
  # /* read and assign file from log list */
  print(beast_phylo)
  assign(paste("log", i, sep=''), read.table(beast_phylo, header=T))
  attach(get(paste("log", i, sep='')))
  
  R0_names = names(get(paste("log", i, sep='')))[which(regexpr("reproductiveNumber.", names(get(paste("log", i, sep=''))))>0)]
  delta_names = names(get(paste("log", i, sep='')))[which(regexpr("becomeUninfectiousRate", names(get(paste("log", i, sep=''))))>0)]
  sampling_names = names(get(paste("log", i, sep='')))[which(regexpr("samplingProportion", names(get(paste("log", i, sep=''))))>0)]
  rhosampling=0
  if (length(sampling_names)==0){
    sampling_names = names(get(paste("log", i, sep='')))[which(regexpr("rho", names(get(paste("log", i, sep=''))))>0)]
    rhosampling=1
  }
  
  treeheights = get(paste("log", i, sep=''))[,match("treeheight", tolower(names(get(paste("log", i, sep='')))))]
  origins = get(paste("log", i, sep=''))$origin
  width=median(origins)
  
  nsamples= length(get(R0_names[1]))
  burnin = round(burninpercent*nsamples/100)
  
  intervalNumber = length(R0_names)
  if (intervalNumber > gridSize) {gridSize = intervalNumber}
  
  medians = matrix(data=NA, nrow= 1, ncol = gridSize)
  medians_G = matrix(data=NA, nrow= 1, ncol = gridSize)
  medians_H = matrix(data=NA, nrow= 1, ncol = gridSize)
  
  hpd_F = matrix(data=NA, nrow= 2, ncol = gridSize)
  hpd_G = matrix(data=NA, nrow= 2, ncol = gridSize)
  hpd_H = matrix(data=NA, nrow= 2, ncol = gridSize)
  
  F = matrix(data=NA, nrow=nsamples-burnin, ncol = gridSize) #reproductiveNumber
  G = matrix(data=NA, nrow=nsamples-burnin, ncol = gridSize) #becomeuninfectiousRate
  H = matrix(data=NA, nrow=nsamples-burnin, ncol = gridSize) #samplingProportion
  
  step = width/(gridSize-1)
  recent = yday(ymd(recent)) /365
  F_times = seq(recent-width, recent, step)
  
  for(k in 1:(nsamples-burnin)){
    # print(k)
    time = origins[k+burnin]
    idx_all  = c()
    for (l in 1:length(F_times)){
      currentWidth = time / intervalNumber
      index = ceiling(intervalNumber - (recent - F_times[l])/currentWidth )
      #	print(index)
      F[k,l] = get(R0_names[max(index,1)])[k+burnin]
      
      #	if (length(deltal_names)==length(R0_names))	
      #						G[k,l] = get(delta_names[max(index,1)])[k+burnin]
      #					else	G[k,l] = get(delta_names[1])[k+burnin]
      if (length(sampling_names)==length(R0_names))	
        H[k,l] = get(sampling_names[max(index,1)])[k+burnin]
      else	H[k,l] = get(sampling_names[1])[k+burnin]
      
    }
  }
  
  for(j in 1:gridSize){
    if (length(which(F[,j]!="NA")) > (nsamples / 10)) {
      medians[1,j] = median(F[,j],na.rm=T)
      #	medians_G[1,j] = median(G[,j],na.rm=T)
      medians_H[1,j] = median(H[,j],na.rm=T)
      hpd_F[,j] = boa.hpd(F[which(F[,j]!="NA"),j], 0.05)[1:2]
      #	hpd_G[,j] = boa.hpd(G[which(G[,j]!="NA"),j], 0.05)[1:2]
      hpd_H[,j] = boa.hpd(H[which(H[,j]!="NA"),j], 0.05)[1:2]
    }
  }
  
  #layout20x(3,1)
  
  # /* plot reproductiveNumber */
  #plot(1, ylab='Reproductive number', xlim=c(recent-width, recent), ylim=c(0,max(hpd_F[2,],na.rm=T)*1.1), xlab="Year", col='white', main = '')
  #m#inor.tick(nx=5, ny=2, tick.ratio=.2)
  #polygon(c(F_times,rev(F_times)), c(hpd_F[2,], rev(hpd_F[1,])),col = "grey90", border = NA)
  #lines(c(F_times), c(medians[1,]), type='l')
  #abline(1,0,col='grey')
  
  
  day = round((F_times* 365))
  dates = as.Date(day,origin="2020-01-01")
  return(data.frame(dates=dates,ro=medians[1,],ro_lo=hpd_F[1,],ro_hi=hpd_F[2,]))
  # p2 = data.frame(dates=dates,ro=medians[1,],ro_lo=hpd_F[1,],ro_hi=hpd_F[2,]) %>% ggplot(aes(y=ro,x=(dates),group=1))+ geom_ribbon(aes(ymin=ro_lo,ymax=ro_hi),fill="grey79")  + geom_point() + geom_line()  + theme_bw()  + 
  #   xlab("Date") +  ylab("Reproductive Number (R0)") + geom_hline(yintercept = 1) + geom_vline(xintercept = ymd("2020-03-20")) +geom_vline(xintercept = ymd("2020-04-07")) +
  #   geom_vline(xintercept = ymd("2020-05-18")) +scale_x_date(date_breaks = "15 days",limits = c(ymd("2020-02-20"),ymd(recent)))
  
  
  
  
  # p1 = gisaid_phased_filt %>% filter(location == "Los Angeles County") %>% ggplot(aes(x=date_fix)) + geom_histogram(binwidth = 7) + scale_x_date(date_breaks = "15 days",limits = c(ymd("2020-02-01"),ymd(recent))) + theme_bw()
  
  
  
  
  
}

all_b_1_43_lockdown = process_beast_phylo("data//b_1_43_lockdown.log","2020-04-20")
all_b_1_43_lockdown %>% ggplot(aes(y=ro,x=(dates),group=1))+ geom_ribbon(aes(ymin=ro_lo,ymax=ro_hi),fill="grey79")  + geom_line()  + theme_bw()  + 
  xlab("Date") +  ylab("Reproductive Number (Re)") + geom_hline(yintercept = 1) + geom_vline(xintercept = ymd("2020-03-20"),color="#e41a1c") +geom_vline(xintercept = ymd("2020-04-07"),color="#377eb8") +
  geom_vline(xintercept = ymd("2020-05-18"),color="#4daf4a") +scale_x_date(date_breaks = "2 weeks") + theme(text=element_text(size=24)) 
ggsave(file="paper/figures/3.png",width=16,height=12,dpi=300)

#all_b_1_43_lockdown = process_beast_phylo("data/beast_phylo/1/b_1_43_lockdown.log","2020-04-20")
#all_b_1_43_lockdown %>% ggplot(aes(y=ro,x=(dates),group=1))+ geom_ribbon(aes(ymin=ro_lo,ymax=ro_hi),fill="grey79")  + geom_line()  + theme_bw()  + 
#  xlab("Date") +  ylab("Reproductive Number (Re)") + geom_hline(yintercept = 1) + geom_vline(xintercept = ymd("2020-03-20"),color="#e41a1c") +geom_vline(xintercept = ymd("2020-04-07"),color="#377eb8") +
#  geom_vline(xintercept = ymd("2020-05-18"),color="#4daf4a") +scale_x_date(date_breaks = "2 weeks") + theme(text=element_text(size=24)) 


#all_la_df3 = process_beast_phylo("data/beast_phylo/la3.log","2020-06-28")
#all_la_lock = process_beast_phylo("data/beast_phylo/focus_on_lockdown.log","2020-04-20")
