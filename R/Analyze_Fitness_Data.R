Analyze_Fitness_Data<- function(Well_key, FC_data){
  
  options(scipen = 999)
  #user supplies boolean values for grouping by replicates (in which case identical 
  #competions in the well key will be grouped), calculating error, or plotting results. 
  group_replicates<- as.logical(readline(prompt="Group replicates? (TRUE/FALSE): "))
  error_method<- as.logical(readline(prompt="Calculate error? (TRUE/FALSE): "))
  plot_boolean<- as.logical(readline(prompt="Plot results? (TRUE/FALSE): "))
  
  #remove rows for Mean and SD from FlowJo table output 
  if("Std.Error" %in% colnames(FC_data)){
    FC_data<-FC_data[-which(colnames(FC_data)=="Std.Error")]
  }
  if ("Mean" %in% FC_data[,1]){
    FC_data<-FC_data[-(which(FC_data[,1]=="Mean")), ]
  }
  if ("SD" %in% FC_data[,1]){
    FC_data<-FC_data[-(which(FC_data[,1]=="SD")), ]
  }
  if (identical(Well_key[,1], FC_data[,1])){
    result_df<- data.frame(Well_key)
    result_df$"Competition"<-paste(Well_key[,2], Well_key[,3])
    result_df<-result_df[c(1,5)]
    result_df$selection_coefficient<-NA
    if(error_method){
      result_df$"Std Err of Regression"<-NA
      result_df$"CI"<-NA
    }
    
    time_points<- (c(0:((ncol(FC_data)-3)/2)))*10
    
    if(!group_replicates){
      for (i in 1:nrow(FC_data)){
        ###grouping identical competitions add here 
        if (Well_key[i,4]==Well_key[i,2]){
          #if data is gated on experimental strain
          exp_vector<- FC_data[i, (seq(3, ncol(FC_data), 2))]
          count_vector<- FC_data[i, (seq(2, ncol(FC_data), 2))]
          ref_vector<- count_vector - exp_vector
          ln_exp_ref_vector<- log(exp_vector/ref_vector)
          ln_exp_ref_vector<-as.numeric(ln_exp_ref_vector[1,])
          linmod <- lm(na.exclude(ln_exp_ref_vector ~ time_points))
          result_df[i,3]<-signif((as.numeric(coef(linmod)[2])), 4)
          if(error_method){
            if (is.na(result_df[i,3])){
              result_df[i,4]<-NA
              result_df[i,5]<-NA
            }
            else{
            result_df[i,4]<-coef(summary(linmod))[2,2]
            crit_t<- qt(1-.05/2,(length(time_points-2)))
            result_df[i,5]<- signif(crit_t*result_df[i,4], 4)
          }
          }
        }
        else { #assumes data is gated on reference
          
          ##grouping identical competitions add here 
          ref_vector<- FC_data[i, (seq(3, ncol(FC_data), 2))]
          count_vector<- FC_data[i, (seq(2, ncol(FC_data), 2))]
          exp_vector<- count_vector - ref_vector
          ln_exp_ref_vector<- log(exp_vector/ref_vector)
          ln_exp_ref_vector<-as.numeric(ln_exp_ref_vector[1,])
          linmod <- lm(na.exclude(ln_exp_ref_vector ~ time_points))
          result_df[i,3]<-signif((as.numeric(coef(linmod)[2])), 4)
        
          if(error_method){
            if (is.na(result_df[i,3])){
              result_df[i,4]<-NA
              result_df[i,5]<-NA
            }
            else{
              result_df[i,4]<-coef(summary(linmod))[2,2]
              crit_t<- qt(1-.05/2,(length(time_points-2)))
              result_df[i,5]<- signif(crit_t*result_df[i,4], 4)
            }
          }
        }
      }
      if (plot_boolean){
        library(ggplot2)
        upper_y<- max(result_df[,3])+ result_df[which.max(result_df[,3]), 5] + .02
        lower_y<-min(result_df[,3]) - result_df[which.min(result_df[,3]), 5] - .02
        if (error_method){
          a<- ggplot(result_df, aes(x=Well.ID, y=selection_coefficient))+
            geom_errorbar(aes(ymin=selection_coefficient-CI, ymax=selection_coefficient+CI), width=.1)+
            geom_point(aes(colour=Competition))+
            ylim(lower_y, upper_y)+
            xlab("Competition")+
            ylab("Fitness advantage")+
            theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
          print(a)
        }
        else {
          a<- ggplot(result_df, aes(x=Competition, y=selection_coefficient))+
            geom_point(aes(colour=Competition))+
            ylim(lower_y, upper_y)+
            xlab("Competition")+
            ylab("Fitness advantage")+
            theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
          print(a)
        }
      }
    }
    if (group_replicates){
      library(reshape2)
      result_df$"replicates"<-NA
      groups<- c(1)
      Well_key<- Well_key[order(Well_key[,2], Well_key[,3]), ]
      FC_data<-FC_data[order(match(FC_data[,1], Well_key[,1])), ]
      result_df<- result_df[order(match(result_df[,1], Well_key[,1])), ]
      for (i in 2:nrow(Well_key)){
        ids<- as.list(Well_key[i, 2:4])
        prev_ids<-as.list(Well_key[i-1, 2:4])
        if (!identical(ids, prev_ids)){
          groups<- append(groups, i)
        }
      }
      for (j in 1:length(groups)){
        if(length(groups)>j){
          end<- groups[j+1]-1
        } else {
          end<-nrow(FC_data)
        }
       
        replicate_group<- FC_data[groups[j]:end, 2:ncol(FC_data)]
        count_df<- replicate_group[,(seq(1, ncol(replicate_group), 2))]
        if (Well_key[groups[j],4]==Well_key[groups[j],2]){ #if gated on experimental
          exp_df<- replicate_group[, (seq(2, ncol(FC_data), 2))]
          ref_df<- count_df - exp_df
        }  else { #assumes gated on reference 
          ref_df<- replicate_group[, (seq(2, ncol(FC_data), 2))]
          exp_df<- count_df - ref_df
        } 
        ln_exp_ref_df<- log(exp_df/ref_df)
        colnames(ln_exp_ref_df)<-time_points
        ln_exp_ref_df<-melt(ln_exp_ref_df)
        ln_exp_ref_df[,1]<-as.numeric(ln_exp_ref_df[,1])
        linmod<-lm(formula= value ~ variable, data=ln_exp_ref_df)
        result_df[groups[j]:end,3]<-signif((as.numeric(coef(linmod)[2])), 4)
        if(error_method){
          result_df[groups[j]:end,4]<-coef(summary(linmod))[2,2]
          crit_t<- qt(1-.05/2,((length(time_points)*(nrow(replicate_group)))-2))
          result_df[groups[j]:end,5]<- signif(crit_t*result_df[groups[j],4], 4)
        }    
        result_df$replicates[groups[j]]<-nrow(replicate_group)
      }
      if(plot_boolean){
        library(ggplot2)
        library(ggrepel)
        upper_y<- max(result_df[,3])+ result_df[which.max(result_df[,3]), 5] + .02
        lower_y<-min(result_df[,3]) - result_df[which.min(result_df[,3]), 5] - .02
        
        if (error_method){
          a<- ggplot(result_df, aes(Competition, selection_coefficient))+
            geom_errorbar(aes(ymin=selection_coefficient-CI, ymax=selection_coefficient+CI), width=.1)+
            geom_point(aes(colour=Competition), size=2)+
            ylim(lower_y, upper_y)+
            xlab("Competition")+
            ylab("Fitness advantage")+
            theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
            geom_text_repel(aes(label=replicates)) 
          print(a) 
        }
        else {
          a<- ggplot(result_df, aes(Competition, selection_coefficient))+
            geom_point(aes(col=Competition), size=2)+
            ylim(lower_y, upper_y)+
            xlab("Competition")+
            ylab("Fitness advantage")+
            theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
            geom_text_repel(aes(label=replicates)) 
          print(a) 
        }
      }
    }
    return (result_df)
  }
  else {
    print("Error - column 1 of well_key and column 1 of data must match")
  }
  
}