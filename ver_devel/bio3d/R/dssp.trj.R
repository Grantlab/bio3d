dssp.trj <- function(pdb, trj, skip=1000, threshold=3, file.head="") {
###Filtering
filter.raw <- seq(1,dim(trj)[1],skip)
####DSSP calculation
dssp.ref<-dssp(pdb)
helix.ref<-sum(dssp.ref$helix$length)/sum(pdb$calpha)*100
sheet.ref<-sum(dssp.ref$sheet$length)/sum(pdb$calpha)*100

unfold.frames <- NULL
for (frame in filter.raw){

	pdb.temp<-pdb
	pdb.temp$xyz<-trj[frame,]
	dssp.test<-dssp(pdb.temp)
	helix.test<-sum(dssp.test$helix$length)/sum(pdb$calpha)*100
	sheet.test<-sum(dssp.test$sheet$length)/sum(pdb$calpha)*100
	helix.diff<-helix.ref - helix.test
	sheet.diff<-sheet.ref - sheet.test
	
	if (helix.diff >= helix.ref/threshold || sheet.diff >= sheet.ref/threshold){

		print (paste('Warning! Check frame', frame, 'in your trajectory! Possible unfolding events, frame saved'))
		write.pdb(xyz=trj[frame,], file=paste(file.head, frame, '.pdb', sep=""),pdb=pdb)
		unfold.frames <- c(unfold.frames, frame)
	} 
}
if(is.null(unfold.frames)) {
	print ('Everything looks good!')
}	
return(unfold.frames)
}
