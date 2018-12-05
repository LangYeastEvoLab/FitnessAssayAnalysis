


Fitness_ANCOVA_post_hoc<- function(ANCOVA_results){
  post_hoc_comparisons<- TukeyHSD(ANCOVA_results, which='Competition')
  a<-plot(post_hoc_comparisons)
  print(a)
  results<-as.data.frame(post_hoc_comparisons$Competition)
  for (i in 1:nrow(results)){
    if (results[i,4]<.01){
      print(paste(rownames(results)[i], "is significant at p = ", results[i,4]))
    }
  }
  return(results)
}