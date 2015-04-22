###########################
##-- BLAST and PDB INPUT  #
###########################

##- PDB input UI that is responsive to 'reset button' below
output$resetable_pdb_input <- renderUI({
  ## 'input$reset_pdb_input' is just used as a trigger for reset
  reset <- input$reset_pdb_input
  textInput("pdbid", label="Enter RCSB PDB code/ID:", value = "4Q21") #)
})


### Input sequence ###
get_sequence <- reactive({

  ## option 1 - PDB code provided
  if(input$input_type == "pdb") {
    anno <<- input_pdb_annotation()
  
    if(is.vector(input$chainId)) {
      ind <- which(anno$chainId==input$chainId[1])
      seq <- unlist(strsplit(anno$sequence[ind], ""))
    }
    else {
      seq <- unlist(strsplit(anno$sequence[1], ""))
    }
    return(as.fasta(seq))
  }

  ## option 2 - sequence provided
  if(input$input_type == "sequence") {
    go <- input$action_input

    inp <- unlist(strsplit(input$sequence, "\n"))
    inds <- grep("^>", inp, invert=TRUE)

    if(!length(inds)>0)
      stop("Error reading input sequence")

    inp <- toupper(paste(inp[inds], collapse=""))

    print(inp)
    
    seq <- as.fasta(unlist(strsplit(inp, "")))
    print(seq)
    return(seq)
  }
})

### Annotate input PDB
input_pdb_annotation <- reactive({
  
  if(nchar(input$pdbid)==4) {
    progress <- shiny::Progress$new(session, min=1, max=5)
    on.exit(progress$close())
    
    progress$set(message = 'Fetching PDB data',
                 detail = 'This should be quick...')
    progress$set(value = 3)

    anno <- get_annotation(input$pdbid, use_chain=FALSE)
    for(i in 4:5)
      progress$set(value = i)
    
    return(anno)
  }
  else {
    stop("Provide a 4 character PDB code")
  }
})

### Run BLAST
run_blast <- reactive({
  message("blasting... ")
  input_sequence <<- get_sequence()
  
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
      hmm <- phmmer_local(input_sequence)
    else 
      hmm <- hmmer(input_sequence)
    
    progress$set(value = 5)
  }
  
  hmmer_data <<- hmm
  return(hmmer_data)
})

filter_hits <- reactive({
  blast <- run_blast()
  hit_inds <- which(filter.hmmer(blast, cutoff=input$cutoff)$inds)

  if(is.null(input$limit_hits)) {
    limit <- 5
  }
  else {
    limit <- as.numeric(input$limit_hits)
  }
  
  if(limit > 0)
    hit_inds <- hit_inds[1:limit]

  n <- length(blast$acc)
  m <- n - length(hit_inds)
  if (m > 5)
    m <- 5

  not_hit_inds <- NULL
  if(m > 0) {
     not_hit_inds <- (hit_inds[length(hit_inds)]+1):(length(hit_inds)+m)
  }

  return(list(hit_inds=hit_inds, not_hit_inds=not_hit_inds, inds=c(hit_inds, not_hit_inds)))
})

get_hit_ids <- reactive({
  blast <- run_blast()
  inds <- filter_hits()
  hits <- blast[inds$hit_inds, "acc"]
  return(hits)
  })

### Blast plot on front page
# output$blast_plot <- renderPlot({
#   blast <- run_blast()
#   hits <- filter.hmmer(blast, cutoff=input$cutoff)
#   plot.hmmer2(blast, hits)
# })
  output$blast_plot <- renderChart2({
    blast <- run_blast()
    blast$mlog.evalue = blast$score
    hits <- filter.hmmer(blast, cutoff=input$cutoff)
    ## Removed call to plot.hmmer2
    ## and moved plotting withing renderChart
    #plot.hmmer2(blast, hits)
    
    if(is.null(hits)) {
      hits <- filter.hmmer(blast)
    }
    cutoff <- hits$cutoff
    gp <- hits$gp.hits
    gps <- hits$gps
    z <- blast$score
    ## generate a dataframe: pc$z + pdbids + group
    data <- data.frame(
      Bitscore = z,
      id = blast$acc,
      group = sapply(gps, function(x) if(x==1) "Above" else "Below cutoff"),
      Hits = 1:length(z)
    )
    
    p1 <- nPlot(
      x = "Hits", y = "Bitscore",
      group = "group",
      data = data,
      type="scatterChart"
    )
    p1$chart( color = c("red","black") )
    p1$xAxis( axisLabel = "Hits" )
    p1$yAxis( axisLabel = "Bitscore", width=50 )
    p1$chart( forceY = if( !(min(data$Bitscore)<0) && min(data$Bitscore)-20 < 0  ) 
      c(0,max(data$Bitscore)+20) else range(data$Bitscore) + c(-20, 20) )
    p1$chart(tooltipContent = "#! function(key, x, y, e){ 
      return '<b>pdb:</b> ' + e.point.id
    } !#")
    p1$setTemplate(afterScript = '<script>
      var css = document.createElement("style");
      css.type = "text/css";
      css.innerHTML = ".nv-y .nv-axislabel { font-size: 30px; }";
      document.body.appendChild(css);
    </script>'
  )
    p1
  
  })

### Table of annotated BLAST results
output$blast_table <- renderDataTable({
  blast <- run_blast()
  inds <- filter_hits()

  hits <- blast[inds$inds,, drop=FALSE]
  anno <- get_annotation(hits$acc)
  
  anno$url <- paste0("<a href=\"", "http://pdb.org/pdb/explore/explore.do?structureId=", substr(anno$acc, 1, 4), "\" target=\"_blank\">", anno$acc, "</a>")
  anno$score <- hits$score
  anno$id <- 1:nrow(anno)

  checked <- rep("CHECKED", length(inds$inds))
  checked[inds$not_hit_inds] <- ""
  checkbox <- paste0("<input type=\"checkbox\" name=\"pdb_ids\" value=\"", hits$acc, "\"",  checked, ">")
  anno$check <- checkbox
  
  return(anno[, c("id", "check", "url", "compound", "source", "ligandId", "chainLength", "score")])
}, escape=FALSE, options = list(lengthChange=FALSE, paging=FALSE))
                                      
  #,                                    
  #                                      callback = "function(table) {
  #    table.on('click.dt', 'tr', function() {
  #      $(this).toggleClass('selected');
  #      Shiny.onInputChange('rows',
  #                          table.rows('.selected').indexes().toArray());
  #    });
  #  }")
                                      



## checkbox 
output$pdb_chains <- renderUI({
  anno <- input_pdb_annotation()
  radioButtons("chainId", label="Choose chain ID:",
               choices=anno$chainId, inline=TRUE)
  
})


output$resetable_cutoff_slider <- renderUI({
  reset <- input$reset_cutoff

  blast <- run_blast()
  cutoff <- filter.hmmer(blast, cutoff=NULL)$cutoff

  sliderInput("cutoff", "Adjust cutoff:",
              min = floor(min(blast$score)), max = floor(max(blast$score)), value = cutoff)
})



#output$cutoff_slider <- renderUI({
#  blast <- blast()##

#  if(is.null(input$cutoff))
#    cutoff <- filter.hmmer(blast, cutoff=NULL)$cutoff
#  else
#    cutoff <- input$cutoff
#  
#  sliderInput("cutoff", "Set cutoff:",
#              min = floor(min(blast$score)), max = floor(max(blast$score)), value = cutoff)
#})


output$hits_slider <- renderUI({
  blast <- run_blast()
  hits <- filter.hmmer(blast, cutoff=input$cutoff, value=TRUE)$acc

  sliderInput("limit_hits", "Limit hits:",
              min = 1, max = length(hits), value = 5, step=1)
})


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

getstuff <- function() {
  blast <- run_blast()
  inds <- filter_hits()

  hits <- blast[inds$inds, c("acc", "score"), drop=FALSE]
  hits$id <- 1:nrow(hits)
  return(hits)
}

blast_plot2 <- renderChart2({
  invisible(capture.output( hits <- getstuff() ))
  m1 <- dPlot(x="id", y="score", z="score", groups="acc", type="bubble", data=hits)
  return(m1)
})

blast_plot3 <- renderChart2({
  d8 <- dPlot(
    x = c("indexname","date"),
    y = "value",
    z = "value",
    groups = "indexname",
    data = ust.melt,
    type = "bubble"
    )
  d8$xAxis( grouporderRule = "date", orderRule = "maturity" )
  d8$zAxis( type = "addMeasureAxis", overrideMax = 10 )
  d8
})
