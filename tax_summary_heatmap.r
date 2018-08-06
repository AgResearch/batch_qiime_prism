library(Heatplus)
library(RColorBrewer)
library("gplots")

get_command_args <- function() {
   args=(commandArgs(TRUE))
   if(length(args)!=3 ){
      #quit with error message if wrong number of args supplied
      print('Usage example : Rscript --vanilla  tax_summary_heatmap.r num_profiles=20 moniker=R1_combined_rep_set_otu_table_L7 datafolder=/dataset/public_invermay_scratch/scratch/UW_Calf_microbiome/qiime_bowtie2_analysis/R1uclust')
      print('args received were : ')
      for (e in args) {
         print(e)
      }
      q()
   }else{
      print("Using...")
      # seperate and parse command-line args
      for (e in args) {
         print(e)
         ta <- strsplit(e,"=",fixed=TRUE)
         switch(ta[[1]][1],
            "datafolder" = datafolder <<- ta[[1]][2],
            "moniker" = moniker <<- ta[[1]][2],
            "num_profiles" = num_profiles <<- as.integer(ta[[1]][2])
         )
      }
   }
}


get_data <- function(moniker) {
   data<-read.table(paste(moniker, ".txt",sep=""), header=TRUE, row.names=1, sep="\t")
   return(data)
}


get_proportions <- function(data) {
   proportions <- data[,FALSE]
   for(column_number in sequence(ncol(data))) {
      total_count = sum(data[,column_number])
      proportions <- cbind(proportions , sapply(data[,column_number], "/", total_count))
   }
   colnames(proportions) <- colnames(data)
   return(proportions)
}



get_self_information <- function(data) {
   self_information <- data[,FALSE]
   for(column_number in sequence(ncol(data))) {
      total_count = sum(data[,column_number])
      approx_zero = min(subset(data, data[,column_number] > 0, select = column_number))/2.0
      self_information <- cbind(self_information,  -log(sapply(data[,column_number],"max",approx_zero)/total_count, 2))
   }
   colnames(self_information) <- colnames(data)
   return(self_information)
}


draw_heatmap <- function(datamatrix, num_clust) {
   # want to plot num_clust  broad taxonomy profiles - so cluster the profiles (if necessary)
   if(nrow(datamatrix) > 1.5 * num_clust) {
      clustering <<- kmeans(datamatrix, num_clust, iter.max=500)


      # label each profile with the name of the species whose profile
      # is closest to the center of each cluster - so find these

      closest_dists = rep(NA,nrow(clustering$centers))
      closest_rownums = rep(NA,nrow(clustering$centers))

      for (center_num in sequence(nrow(clustering$centers))) {
         v_center = as.numeric(clustering$centers[center_num,])
         for (row_num in sequence(nrow(datamatrix))) {
            # only consider this row as potentially supplying a name for the cluster if it is in the cluster
            # (sometimes a point not in the cluster can be closer to the center than a point that is in the cluster)
            if ( clustering$cluster[row_num] == center_num ) {
               v_data = as.numeric(datamatrix[row_num,])

               # calculate the distance from the center and update the closest_dists data structure
               d = (v_center - v_data) %*% (v_center - v_data)
               if(is.na(closest_dists[center_num])) {
                  closest_dists[center_num] = d
                  closest_rownums[center_num] = row_num
               }
               else if( d < closest_dists[center_num] ) {
                  closest_dists[center_num] = d
                  closest_rownums[center_num] = row_num
               }
            } 
         }
      }

      # assign the labels to the clustered data
      rownames=rownames(datamatrix)[closest_rownums]
      clustered_data = clustering$centers
      rownames(clustered_data) = rownames
   }
   else {
      clustered_data = datamatrix
      clustering <<- NA
   }

   # ref for configuring plot
   #http://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
   #1 Heatmap,
   #2 Row dendrogram,
   #3 Column dendrogram,
   #4 Key


   # draw the heatmap in the usual way
   #cm<-brewer.pal(11,"Spectral") # a diverging palette
   cm<-brewer.pal(9,"OrRd") # a sequential palette 
   cm <- rev(cm)


   # set up a vector which will index the labels that are to be blanked out so that 
   # only every nth col is labelled, 
   # the rest empty strings, n=col_label_interval.
   number_of_column_labels=40
   col_label_interval=max(1, floor(ncol(clustered_data)/number_of_column_labels))  # 1=label every location 2=label every 2nd location  etc 
   colLabels <- colnames(as.matrix(clustered_data))
   colBlankSelector <- sequence(length(colLabels))
   colBlankSelector <- subset(colBlankSelector, colBlankSelector %% col_label_interval != 0) 
                       # e.g. will get (2,3, 5,6, 8,9, ..)
                       # so we will only label rows 1,4,7,10,13 etc)


   # initial plot to get the column re-ordering
   jpeg(filename = paste(moniker, ".jpg",sep="") , width=830, height=1200) # with dendrograms

   hm<-heatmap.2(as.matrix(clustered_data),  scale = "none", 
   #hm<-heatmap.2(as.matrix(datamatrix),  scale = "none", 
       dendrogram = "col",  
       trace="none",
       #trace = "none", breaks =  -2 + 4/9*seq(0,11), 
       col = cm , key=FALSE, density.info="none", 
       #keysize=1.0, margin=c(17,25), cexRow=1.5, cexCol=1.6, 
       keysize=1.0, margin=c(17,28), cexRow=1.5, cexCol=1.6, 
       lmat=rbind(  c(4,3,0 ), c(2, 1, 0) ), lwid=c(.2, .6, 0 ), lhei=c(.5, 3))

  dev.off()


   # edit the re-ordered vector of col labels, obtained from the heatmap object, so that only 
   # every nth label on the final plot has a non-empty string
   # this is for the internal distance matrix
   indexSelector <- hm$colInd[length(hm$colInd):1]    
   indexSelector <- indexSelector[colBlankSelector]
   colLabels[indexSelector] = rep('',length(indexSelector))

   jpeg(filename = paste(moniker, ".jpg",sep=""), width=1400, height=1000) # with dendrograms
   hm<-heatmap.2(as.matrix(clustered_data),  scale = "none", 
       dendrogram = "col",  
       trace="none",
       #trace = "none", breaks =  -2 + 4/9*seq(0,11), 
       col = cm , key=FALSE, density.info="none", 
       #keysize=1.0, margin=c(17,25), cexRow=1.5, cexCol=1.6, 
       #keysize=1.0, margin=c(27,28), cexRow=1.2, cexCol=1.2, 
       #keysize=1.0, margin=c(27,48), cexRow=1.2, cexCol=1.2, 
       keysize=1.0, margin=c(27,78), cexRow=1.3, cexCol=1.3, 
       #lmat=rbind(  c(4,3,0 ), c(2, 1, 0) ), lwid=c(.2, .6, 0 ), lhei=c(.25, 3),labCol=colLabels)
       lmat=rbind(  c(4,3,0), c(2,1,0)), lwid=c(.1, 1.2, 0), lhei=c(.25, 3 ),labCol=colLabels)

   # the column labels on the plots are usually too crowded so supply a file with the 
   # column names ordered as per the plot
   write.table(colnames(as.matrix(clustered_data))[hm$colInd[1:length(hm$colInd)]] , file=paste(moniker, "_samplenames_ordered.dat",sep=""),row.names=TRUE,sep="\t")
   # the row labels on the plots may be truncated  so supply a file with the 
   # row  names ordered as per the plot
   write.table(rownames(as.matrix(clustered_data))[hm$rowInd[length(hm$rowInd):1]] , file=paste(moniker, "_taxnames_ordered.dat",sep=""),row.names=TRUE,sep="\t")

   if ( ! is.na( clustering ) ) {
      # supply the tax clusters
      write.table(clustering$cluster, file=paste(moniker, "_tax_clusters.dat",sep=""),row.names=TRUE,sep="\t")
      # supply the names given to the tax clusters
      write.table(rownames, file=paste(moniker, "_tax_cluster_names.dat",sep=""),row.names=TRUE,sep="\t")
   }


   # 
   clust = as.hclust(hm$colDendrogram)
   sink(paste(moniker, "_heatmap_clustering_support.txt",sep=""))
   print("clust$merge:")
   print(clust$merge)
   print("clust$height:")
   print(clust$height)
   print("clust$order")
   print(clust$order)
   print("clust$labels")
   print(clust$labels)
   sink()
   write.table(cutree(clust, 1:dim(clustered_data)[2]),file=paste(moniker, "_heatmap_clusters.txt",sep=""),row.names=TRUE,sep="\t")  # ref https://stackoverflow.com/questions/18354501/how-to-get-member-of-clusters-from-rs-hclust-heatmap-2


   dev.off()
}


draw_abundant_heatmap <- function(datamatrix) {

   # want to plot the 40 most abdundant taxa (no clustering).
   # order the data by the total of each row (append the row totals
   # as a column and sort on that)
   sdatamatrix <- cbind(datamatrix, rowSums(datamatrix))
   sdatamatrix <- sdatamatrix[order(sdatamatrix[,ncol(sdatamatrix)]),]    # nb here smaller means more abundant as these are self-information measures 
   sdatamatrix <- head(sdatamatrix, 40)                    # take the first 40 
   sdatamatrix <- sdatamatrix[, sequence(ncol(sdatamatrix)-1)]   # drop the totals column

   # draw the heatmap in the usual way
   #cm<-brewer.pal(11,"Spectral") # a diverging palette
   cm<-brewer.pal(9,"OrRd") # a sequential palette 
   cm <- rev(cm)

   # set up a vector which will index the labels that are to be blanked out so that 
   # only every nth col is labelled, 
   # the rest empty strings, n=col_label_interval.
   number_of_column_labels=40
   col_label_interval=max(1, floor(ncol(sdatamatrix)/number_of_column_labels))  # 1=label every location 2=label every 2nd location  etc 
   colLabels <- colnames(as.matrix(sdatamatrix))
   colBlankSelector <- sequence(length(colLabels))
   colBlankSelector <- subset(colBlankSelector, colBlankSelector %% col_label_interval != 0) 
                       # e.g. will get (2,3, 5,6, 8,9, ..)
                       # so we will only label rows 1,4,7,10,13 etc)


   # initial plot to get the column re-ordering
   jpeg(filename = "hm_internal.jpg" , width=830, height=1300) # with dendrograms

   hm<-heatmap.2(as.matrix(sdatamatrix),  scale = "none", 
       dendrogram = "col",  
       trace="none",
       #trace = "none", breaks =  -2 + 4/9*seq(0,11), 
       col = cm , key=FALSE, density.info="none", 
       #keysize=1.0, margin=c(17,25), cexRow=1.5, cexCol=1.6, 
       keysize=1.0, margin=c(17,28), cexRow=1.5, cexCol=1.6, 
       lmat=rbind(  c(4,3,0 ), c(2, 1, 0) ), lwid=c(.2, .6, 0 ), lhei=c(.5, 3))

  dev.off()


   # edit the re-ordered vector of col labels, obtained from the heatmap object, so that only 
   # every nth label on the final plot has a non-empty string
   # this is for the internal distance matrix
   indexSelector <- hm$colInd[length(hm$colInd):1]    
   indexSelector <- indexSelector[colBlankSelector]
   colLabels[indexSelector] = rep('',length(indexSelector))

   jpeg(filename = paste(moniker, "_abundant.jpg",sep=""), width=1400, height=1000) # with dendrograms
   hm<-heatmap.2(as.matrix(sdatamatrix),  scale = "none", 
       dendrogram = "col",  
       trace="none",
       #trace = "none", breaks =  -2 + 4/9*seq(0,11), 
       col = cm , key=FALSE, density.info="none", 
       #keysize=1.0, margin=c(17,25), cexRow=1.5, cexCol=1.6, 
       keysize=1.0, margin=c(27,78), cexRow=1.3, cexCol=1.3, 
       lmat=rbind(  c(4,3,0 ), c(2, 1, 0) ), lwid=c(.1, 1.2, 0 ), lhei=c(.25, 3),labCol=colLabels)
   title(main="40 most abundant taxa")

   # the column labels on the plots are usually too crowded so supply a file with the
   # column names ordered as per the plot
   write.table(colnames(as.matrix(sdatamatrix))[hm$colInd[1:length(hm$colInd)]] , file=paste(moniker, "_abundant_samplenames_ordered.dat",sep=""),row.names=TRUE,sep="\t")
   # the row labels on the plots may be truncated  so supply a file with the
   # row  names ordered as per the plot
   write.table(rownames(as.matrix(sdatamatrix))[hm$rowInd[length(hm$rowInd):1]] , file=paste(moniker, "_abundant_taxnames_ordered.dat",sep=""),row.names=TRUE,sep="\t")


   dev.off()


   # 
   clust = as.hclust(hm$colDendrogram)
   sink(paste(moniker, "_abundant_heatmap_clustering_support.txt",sep=""))
   print("clust$merge:")
   print(clust$merge)
   print("clust$height:")
   print(clust$height)
   print("clust$order")
   print(clust$order)
   print("clust$labels")
   print(clust$labels)
   sink()
   write.table(cutree(clust, 1:dim(sdatamatrix)[2]),file=paste(moniker, "_abundant_heatmap_clusters.txt",sep=""),row.names=TRUE,sep="\t")  # ref https://stackoverflow.com/questions/18354501/how-to-get-member-of-clusters-from-rs-hclust-heatmap-2

}


main <- function() {
   get_command_args()
   setwd(datafolder)
   data <- get_data(moniker)
   props <- get_proportions(data)
   self_info <- get_self_information(data)
   write.table(props,file=paste(moniker,"_proportions.dat",sep=""),row.names=TRUE,sep="\t")
   write.table(self_info,file=paste(moniker, "_self_information.dat",sep=""),row.names=TRUE,sep="\t")
   draw_heatmap(self_info, num_profiles)
   draw_abundant_heatmap(self_info)
}


main()


