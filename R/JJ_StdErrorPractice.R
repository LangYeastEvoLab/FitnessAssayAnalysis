get_some_stderror<- function(Well_key, FC_data){

  remove_samp<- as.logical(readline(prompt="Eliminate poor samples? (TRUE/FALSE): "))

  #Katie added some code to remove mean and SD rows from flowjo data and to ensure that the well key matches the FC data
  if ("Mean" %in% FC_data[,1]){
    FC_data<-FC_data[-(which(FC_data[,1]=="Mean")), ]
  }
  if ("SD" %in% FC_data[,1]){
    FC_data<-FC_data[-(which(FC_data[,1]=="SD")), ]
  }
  if (!identical(Well_key[,1], FC_data[,1])){
    print("ERROR - first columns of well key and FC data must match")
  }
  
  FC_data$"Std.Error" <- NA #Creates a new column
  c <- ncol(FC_data)

  time_points<- (c(0:((ncol(FC_data)-4)/2)))*10 #finds the number of time points in dataset 
  
  #Goes through dataframe and calculates Std. Error for each sample
  list_for_plots<- list()
  
  for (i in 1:nrow(FC_data)) {
    count<- FC_data[i, seq(2,ncol(FC_data)-1,2)] #-1 so we don't count Std. Error col
    refcount<- FC_data[i, seq(3,ncol(FC_data)-1,2)]
    diff<- count-refcount
    natlog<- log(diff/refcount)
    natlog<-as.numeric(natlog[1,])
    list_for_plots<- append(list_for_plots, natlog) # storing the natlog data for later
    regress<- lm(na.exclude(natlog ~ time_points))
    slope<- as.numeric(coef(regress)[2])
    stderror<- as.numeric(coef(summary(regress))[2,2])
    FC_data[i,c] <- stderror
  }
  
  ## this section plots raw data and stores it to a pdf file in the directory the script is run in
  min.y<- min(na.omit(as.numeric(list_for_plots)))
  max.y<- max(na.omit(as.numeric(list_for_plots)))
  if (abs(min.y)>abs(max.y)){
    max.y<- abs(min.y)
  } else {
    min.y<- -(max.y)
  }
  samples<- seq(from=1, to=length(list_for_plots), by = length(time_points)) #a list to interate by
  pdf("raw_data_plots.pdf", height = (ncol(FC_data)*7), width = 16)
  par(mfrow=c(nrow(FC_data)/4, 4)) #control the margins of the plots
  par(mar=c(4,4,4,4))
  for (i in 1:length(samples)) {
    start<- samples[i]
    end<- samples[i]+(length(time_points)-1)
    plot(time_points, list_for_plots[start:end], pch= 16, cex= 2, type = "o", 
         col= sample(1:600,1,replace = T), ylim = c(min.y, max.y), xlab=NA, ylab=NA)#generates plots
    legend("topleft", legend=paste("std err =", signif(FC_data[i, ncol(FC_data)], 4)))
    title(main = FC_data[i,1], xlab = "Generation", ylab = "ln(exp/ref)")
  }
  dev.off()
  
  pdf("std_error_histogram.pdf")
  hist(FC_data[,ncol(FC_data)], xlab="standard error", ylab="# samples", main="Histogram of standard errors", col="royalblue")
  dev.off()
  
  # this section asks a user for a standard error cuttoff value 
  # and removes samples that don't meet the cuttoff
  good_samples.frame<- FC_data #copies FC dataframe to new frame 
  for (i in 1:ncol(good_samples.frame)){
    good_samples.frame[,i]=NA
  }
  good_samples.frame<-na.omit(good_samples.frame)
  if (remove_samp){
    print("run adjust_well_key() to update well key for deleted samples")
    cutoff<-as.numeric(readline(prompt="Insert cutoff threshold: "))
    
    for (i in 1:nrow(FC_data)){
      if(FC_data[i,c] < cutoff){
        good_samples.frame <- rbind(good_samples.frame, FC_data[i, ])
        }
    }
    FC_data<- good_samples.frame
  }
  return(FC_data)
}

