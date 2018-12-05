
Adjust_well_key<-function(Well_key, Reduced_FC_data){
  j<- which(!Well_key[,1] %in% Reduced_FC_data[,1])
  Well_key<-Well_key[-j, ]
}