

Fitness_ANCOVA<- function(Well_key, FC_data){
  library(reshape2)
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
  Competition<-paste(Well_key[,2], Well_key[,3])
  
  if (identical(Well_key[,1], FC_data[,1])){ #the rest of the code will only execute if the first column of the well key 
    # matches the first column of the flow data. 

 time_points<- (c(0:((ncol(FC_data)-3)/2)))*10 #finds how many time points are in competition
    
ANCOVA_df<- data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c("Competition", "Generation", "ln_exp_ref"))), stringsAsFactors=FALSE)
  
   groups<- c(1)
   Well_key<- Well_key[order(Well_key[,2], Well_key[,3]), ]
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
        FC_data<-FC_data[order(match(FC_data[,1], Well_key[,1])), ]
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
        ln_exp_ref_df<- cbind(Competition[groups[j]:end], ln_exp_ref_df)
        colnames(ln_exp_ref_df)<-c("Competition", "Generation", "ln_exp_ref")
        ANCOVA_df<- rbind(ANCOVA_df, ln_exp_ref_df)
      }
   ANCOVA_results<- aov(ln_exp_ref~Generation*Competition, data=ANCOVA_df)
   summary(ANCOVA_results)[[1]][5][[1]][3]->p_val
   if(p_val<.01){
    print(paste("There are significant differences between competitions at p = ", p_val))
    print("To identify pairwise differences, run Fitness_ANCOVA_post_hoc")
   }
    return (ANCOVA_results)
  }
  else {
    print("Error - column 1 of well_key and column 1 of data must match")
  }
}
