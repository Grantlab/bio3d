########################
##-- BLAST / HMMER     #
########################

blast <- reactive({
  message("blasting")
  seq <- final_sequence()
  
  if(toupper(input$pdbid)=="4Q21" & input$input_type == "pdb")
    load(file="4Q21_hmmer.RData")
  else {
    progress <- shiny::Progress$new(session, min=1, max=5)
    on.exit(progress$close())
    
    progress$set(message = 'Blasting',
                 detail = 'This can be slow...')
    progress$set(value = 2)

    ## use local version if possible
    if(configuration$hmmer$local)
      hmm <- phmmer_local(seq)
    else 
      hmm <- hmmer(seq)
    
    progress$set(value = 5)
  }
  
  hmm <<- hmm
  return(hmm)
})

output$blast_plot <- renderPlot({
  blast <- blast()
  hits <- filter.hmmer(blast, cutoff=input$cutoff)
  plot.hmmer2(blast, hits)
})

output$blast_table <- renderDataTable({
  blast <- blast()
  hits <- blast[filter.hmmer(blast, cutoff=input$cutoff)$inds, , drop=FALSE]

  if(!is.null(input$limit_hits)) {
    n <- as.numeric(input$limit_hits)
    print(n)
    if(n > 0)
      hits <- hits[1:n, , drop=FALSE]
  }
  
  anno <- get_annotation(hits$acc)
  
  url <- paste0("<a href=\"", "http://pdb.org/pdb/explore/explore.do?structureId=", substr(anno$acc, 1, 4), "\" target=\"_blank\">", anno$acc, "</a>")
  anno <- cbind(anno, url)
  anno <- cbind(anno, score=hits$score)
  anno <- cbind(anno, id=1:nrow(anno))

  #checkbox <- paste0("<input type=\"checkbox\" name=\"pdb_ids\" value=\"", hits$acc, "\" CHECKED>")
  #anno <- cbind(anno, check=checkbox)
  
  return(anno[, c("id", "url", "compound", "source", "ligandId", "chainLength", "score")])
},
                                      escape=FALSE
  #,                                    
  #                                      callback = "function(table) {
  #    table.on('click.dt', 'tr', function() {
  #      $(this).toggleClass('selected');
  #      Shiny.onInputChange('rows',
  #                          table.rows('.selected').indexes().toArray());
  #    });
  #  }"
                                      )
                                      

hits1 <- function() {
  blast <- blast()
  hits <- filter.hmmer(blast, cutoff=input$cutoff, value=TRUE)$acc
  
  if(!is.null((input$limit_hits))) {
    n <- as.numeric(input$limit_hits)
    if(n > 0) {
      if(n>length(hits))
        n <- length(hits)
      hits <- hits[1:n]
    }
  }
  return(hits)
}

output$hits2 <- renderPrint({
  invisible(capture.output( hits <- hits1()) )
  print(hits)
})

output$cutoff_slider <- renderUI({
  blast <- blast()

  if(is.null(input$cutoff))
    cutoff <- filter.hmmer(blast, cutoff=NULL)$cutoff
  else
    cutoff <- input$cutoff
  
  sliderInput("cutoff", "Set cutoff:",
              min = floor(min(blast$score)), max = floor(max(blast$score)), value = cutoff)
})


output$hits_slider <- renderUI({
  blast <- blast()
  hits <- filter.hmmer(blast, cutoff=input$cutoff, value=TRUE)$acc

  sliderInput("limit_hits", "Limit hits:",
              min = 1, max = length(hits), value = 5, step=1)
})


output$pdbids_checkboxgroup <- renderUI({
  hits <- hits1()
  names(hits) <- hits
  checkboxGroupInput("pdb_ids", "",
                     hits, selected=hits, inline=TRUE)
})



filter.hmmer <- function(x, cutoff=NULL, cut.seed=NULL, cluster=TRUE, value=FALSE) {
  
  if(is.null(x$evalue))
    stop("missing evalues")

  x$mlog.evalue=x$score

  ##- Find the point pair with largest diff evalue
  dx <- abs(diff(x$mlog.evalue))
  dx.cut = which.max(dx)
  
  
  if(!is.null(cutoff)) {
    ##- Use suplied cutoff
    gps = rep(2, length(x$mlog.evalue))
    gps[ (x$mlog.evalue >= cutoff) ] = 1

  } else {

    if(cluster) {
      ## Ask USER whether to continue with clustering with many hits  
      nhit <- length(x$mlog.evalue)
      if(nhit > 1500) {
        cluster <- readline( paste0(" Note: ", nhit, 
          " hits, continue with TIME-CONSUMING clustering [y/n/q](n): ") )

        cluster <- switch(cluster, y=TRUE, yes=TRUE, q="QUIT", FALSE)
        if(cluster=="QUIT") { stop("user stop") }
      }
    }

    if(is.null(cut.seed)) {
      ## Use mid-point of largest diff pair as seed for
      ##  cluster grps (typical PDB values are ~110)
      cut.seed = mean( x$mlog.evalue[dx.cut:(dx.cut+1)] )
    }

    if(cluster){
      ##- Partition into groups via clustering 
      ##  In future could use changepoint::cpt.var
      hc <- hclust( dist(x$mlog.evalue) )
      if(!is.null(cutoff)) { cut.seed=cutoff } 
      gps <- cutree(hc, h=cut.seed)
    } 

    if(!cluster || (length(unique(gps))==1)) {
      ##- Either we don't want to run hclust or hclust/cutree 
      ##   has returned only one grp so here we will divide   
      ##   into two grps at point of largest diff
      gps = rep(2, length(x$mlog.evalue))
      gps[1:dx.cut]=1
    }
  }

  gp.inds <- na.omit(rle2(gps)$inds)
  gp.nums <- x$mlog.evalue[gp.inds]


  
  cat("  * Possible cutoff values:   ", floor(gp.nums), "\n",
      "           Yielding Nhits:   ", gp.inds, "\n\n")

  if( is.null(cutoff) ) {
    ## Pick a cutoff close to cut.seed
    i <- which.min(abs(gp.nums - cut.seed))
    cutoff <- floor( gp.nums[ i ] )
  }

  inds <- x$mlog.evalue >= cutoff
  cat("  * Chosen cutoff value of:   ", cutoff, "\n",
      "           Yielding Nhits:   ", sum(inds), "\n")
  
  if(value) {
    return(x[inds, ])
  }
  else {
    out <- list(inds=inds, gp.inds=gp.inds, gps=gps, cutoff=cutoff)
    return(out)
  }
}

plot.hmmer2 <- function(x, inds=NULL, cex=1.1) {
  
  if(is.null(inds)) {
    inds <- filter.hmmer(x)
  }

  cutoff <- inds$cutoff
  gp <- inds$gp.inds
  gps <- inds$gps
  z <- x$score
  
  plot(z, xlab="", ylab="Bitscore", col=gps)
  abline(v=gp, col="gray70", lty=3)
  
  pos=c(rep(3, length(gp))[-length(gp)],2)
  text(  gp, z[gp], 
       labels=paste0("Nhit=",gp ,", x=", round(z[gp])), 
       col="black", pos=3, cex=cex)
  
}
