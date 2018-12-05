


carly_row_sums = colSums(carly.dge)
row_sums = merge(bijective, carly_row_sums, by = "row.names")
d <- ggplot(row_sums, aes(x = row_sums["xcoord"], y = row_sums["ycoord"], z = (row_sums["y"])))
d + stat_summary_2d(binwidth = c(156.25,156.25), fun = "sum") + scale_fill_gradientn(colors = c(low = "blue", mid = "yellow",high = "red"))
d + stat_summary_2d(fun = "sum") + scale_fill_gradientn(colors = c(low = "blue", mid = "yellow",high = "red"))

#p = d + stat_summary_2d(binwidth = c(13,13), fun = "sum") + scale_fill_gradientn(colors = c(low = "blue", mid = "yellow",high = "red"))
p = d + stat_summary_2d(binwidth = c(156.25,156.25), fun = "sum") + scale_fill_gradientn(colors = c(low = "blue", mid = "yellow",high = "red"))
p
q=layer_data(p,1)
hist(q["value"])
mn = sapply(q["value"],mean)
s = sapply(q["value"],sd)
legend(x = "topright",legend = c(paste0("sd = ",s) ,paste0("mean = ",mn)))
abline(v = mn, col = "red", lwd = 2)



## This is for making the plot with each of the TSNE assignments plotted in actual space assigning a color to each bead
identity = fpcopy@ident
clustspace = merge(bijective,identity,by="row.names") #bc_locations is generated from the script format_spatial_data
d <- ggplot(clustspace, aes(x = xcoord, y = ycoord, color = y))
d+geom_point() #this will plot the TSNE clusters in space!
#this makes a plot of all of the TSNE clusters where the cluster i will be shown as red and all others as background
# the idendity vector is a vector of all the TSNE idendities with row names as the barcodes. generated above.
for (i in 1:max(as.numeric(identity))){
  print(d+geom_point(data=clustspace, aes(color = ifelse(as.numeric(y)==i, paste0(" ",i-1), paste0("not ",i-1)))))
}

#multiple color plot for the clusters



for (i in all.markers) {
  markers.use = rownames(i)[1:4]
  FeaturePlot(fpcopy,markers.use, reduction.use='spatial')
  
}
