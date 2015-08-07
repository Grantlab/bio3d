###################################
##-- User path for storing stuff  #
###################################
data_path <- reactive({
  dir <- paste0(format(Sys.time(), "%Y-%m-%d"), "_", randstr())
  path <- paste0(configuration$user_data, "/", dir)
  dir.create(path)
  return(path)
})


###########################
##-- PDB AND BLAST INPUT  #
###########################

### set default vaules
rv <- reactiveValues()
rv$pdbid <- "2LUM"
rv$chainid <- "A"
rv$blast <- readRDS("2LUM_blast.RDS")
##rv$pfam <- readRDS("2LUM_pfam.RDS")
rv$limit_hits <- 5
rv$cutoff <- 41
rv$sequence <- "MQYKLVINGKTLKGETTTKAVDAETAEKAFKQYANDNGVDGVWTYDDATKTFTVTE"
rv$pdb_codes <- "1TND, 1KJY_A"

observeEvent(input$pdbid, {
  if(nchar(input$pdbid)>3) {
    if(input$input_type == "pdb") {
      pdbid <- trim(input$pdbid)
      rv$pdbid <- substr(pdbid, 1, 4)
      
      chains <- get_chainids()
      rv$chainid <- chains[1]
      
      if(rv$pdbid == "2LUM") {
        blast <- rv$blast
      }
      else {
        blast <- run_blast()
        rv$blast <- blast
        ##rv$pfam <- pdb.pfam(rv$pdbid)
        ##rv$pfam <- get_pfam_annotation()
      }

      cut <- set_cutoff(blast, cutoff=NULL)
      rv$cutoff <- cut$cutoff
      rv$limit_hits <- 5
      rv$modes <- NULL
      rv$aligned <- FALSE
      rv$fitted <- FALSE
    }
  }
})

observeEvent(input$sequence, {
  if(nchar(input$sequence)>10) {
    if(input$input_type == "sequence") {
      rv$sequence <- input$sequence
      blast <- run_blast()
      rv$blast <- blast
      cut <- set_cutoff(blast, cutoff=NULL)
      rv$cutoff <- cut$cutoff
      rv$limit_hits <- 5
      rv$modes <- NULL
      rv$aligned <- FALSE
      rv$fitted <- FALSE
    }
  }
})

observeEvent(input$chainId, {
  rv$chainid <- input$chainId
})

observeEvent(input$limit_hits, {
  rv$limit_hits <- as.numeric(input$limit_hits)
  rv$modes <- NULL
  rv$aligned <- FALSE
  rv$fitted <- FALSE
})

observeEvent(input$cutoff, {
  rv$cutoff <- as.numeric(input$cutoff)
  rv$modes <- NULL
  rv$aligned <- FALSE
  rv$fitted <- FALSE
})

observeEvent(input$reset_cutoff, {
  if(input$input_type != "multipdb") {
    blast <- rv$blast
    cut <- set_cutoff(blast, cutoff=NULL)

    updateSliderInput(session, "cutoff", value = cut$cutoff)
    updateSliderInput(session, "limit_hits", value = 5)
    
    rv$modes <- NULL
    rv$aligned <- FALSE
    rv$fitted <- FALSE
  }
})

observeEvent(input$reset_pdbid, {
  updateTextInput(session, "pdbid", value = "2LUM")
})


## Returns HMMER results
## dataframe with columns: acc, evalue, score, desc
run_blast <- reactive({
  message("run_blast called")
  ptm <- proc.time()

  
  if(input$input_type != "multipdb") {
    message("blasting")

    if(input$input_type == "pdb")
      if(is.null(rv$chainid))
        return(NULL)

    if(input$input_type == "sequence")
      if(!nchar(rv$sequence) > 10)
        return(NULL)

    input_sequence <- get_sequence()

    if(is.null(input_sequence))
      return(NULL)

    progress <- shiny::Progress$new(session, min=1, max=5)
    on.exit(progress$close())

    progress$set(message = 'Blasting',
                 detail = 'Please wait')
    progress$set(value = 2)

    ## use local version if possible
    if(configuration$hmmer$local)
      hmm <- phmmer_local(input_sequence)
    else
      hmm <- hmmer(input_sequence)

    if(!nrow(hmm) > 0)
      stop("No BLAST hits found")

    ## to uppercase_untouched
    hmm$acc <- format_pdbids(hmm$acc)
    ##blast <- hmm
    ##saveRDS(blast, file="2LUM_blast.RDS")

    progress$set(value = 5)
    progress$close()
    message("time used for blasting")
    print(proc.time() - ptm)
    return(hmm)
  }
  else {
    return()
  }
})


## filters the blast results based on the cutoff
## returns logical vectors of all and limited hits
## as well as accession ids
filter_hits <- reactive({
  message("filter_hits called")

  ##blast <- run_blast()
  blast <- rv$blast
  cutoff <- rv$cutoff

  ## logical vector
  inds <- blast$score > cutoff

  ## limited by input$limit_hits
  limit <- rv$limit_hits
  inds2 <- inds
  if(limit > 0 & limit < sum(inds)) {
    inds2[(limit+1):length(inds2)] <- FALSE
  }

  ## accession ids above cutoff
  hits <- blast$acc[inds]
  hits2 <- blast$acc[inds2]

  out <- list(hits=hits, inds=inds,
              hits_limited=hits2, inds_limited=inds2,
              cutoff=cutoff)

  return(out)
})


set_cutoff <- function(blast, cutoff=NULL) {

  if(is.null(blast))
    return(NULL)

  x <- blast
  cluster <- FALSE
  cut.seed <- NULL

  if(is.null(x$evalue))
    stop("missing evalues")

  x$mlog.evalue <- x$score

  ##- Find the point pair with largest diff evalue
  dx <- abs(diff(x$mlog.evalue))
  dx.cut = which.max(dx)


  if(!is.null(cutoff)) {
    ##- Use suplied cutoff
    gps = rep(2, length(x$mlog.evalue))
    gps[ (x$mlog.evalue >= cutoff) ] = 1

  }
  else {
    if(cluster) {
      ## avoid clustering with many hits
      nhit <- length(x$mlog.evalue)
      if(nhit > 2000) {
        cluster <- FALSE
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

  #cat("  * Possible cutoff values:   ", floor(gp.nums), "\n",
  #    "           Yielding Nhits:   ", gp.inds, "\n\n")

  if( is.null(cutoff) ) {
    ## Pick a cutoff close to cut.seed
    i <- which.min(abs(gp.nums - cut.seed))
    cutoff <- floor( gp.nums[ i ] )
  }

  inds <- x$mlog.evalue >= cutoff
  cat("  * Chosen cutoff value of:   ", cutoff, "\n",
      "           Yielding Nhits:   ", sum(inds), "\n")

  out <- list(inds=inds, gp.inds=gp.inds, grps=gps, cutoff=cutoff)
  return(out)
}



output$blast_plot <- renderUI({
  ##blast <- run_blast()
  blast <- rv$blast
  plotOutput("blast_plot1")

  #if(length(blast$acc) > 100) {
  #  plotOutput("blast_plot1")
  #}
  #else {
  #  showOutput("blast_plot2", "nvd3")

    #tags$script(HTML(
    #  'var css = document.createElement("style");
    #                  css.type = "text/css";
    #                  css.innerHTML = ".nv-x .nv-axislabel { font-size: 20px; }";
    #                  document.body.appendChild(css);
    #                  css = document.createElement("style");
    #                  css.type = "text/css";
    #                  css.innerHTML = ".nv-y .nv-axislabel { font-size: 20px; }";
    #                  document.body.appendChild(css);'
    #  ))
  #}

})

### Blast plot on front page
output$blast_plot1 <- renderPlot({
  ptm <- proc.time()

  
  ##blast <- run_blast()
  blast <- rv$blast
  cut <- set_cutoff(blast, rv$cutoff)

  gp <- cut$gp.inds
  cutoff <- cut$cutoff
  grps <- cut$grps
  z <- blast$score

  if(length(blast$acc) < 650) {
    pch <- 21
    bg <- col
    col <- "grey10"
  }
  else {
    pch <- 1
    bg <- NULL
    col <- col
  }

  col.bg <- sapply(grps, function(x) if(x==1) 'red3' else if(x==2) 'grey50')
  pch = rep(21, length(col.bg))
  pch[col.bg=='grey50'] <- 21
  col.pch = rep('red3', length(pch))
  col.pch[pch==21] <- 'black'


  par(mar=c(6, 4, 0, 0))
  plot(z, xlab="", ylab="Bitscore of Alignment to Input Structure", bg=col.bg, col=col.bg,#col=col.pch, 
    pch=pch, cex=1.1, xaxt='n', cex.axis=1.2)
  axis(1, at=if(length(z) > 20 ) c(1, seq(20, length(z), 20), if(length(z)%%20 > 9) length(z) else NA) else seq(1,20,2), cex.axis=1.2)
  abline(v=gp, col="gray50", lty=3)
  abline(h=cutoff, col="red", lty=2)
  if(!length(input$blast_table_rows_selected)>0) {
    limit <- sort(as.numeric(rv$limit_hits))
    if(limit > 0) {
      points(z[1:limit], col='navyblue', pch=1, cex=2.25)
    }
  }
  else {
    x = sort(as.numeric(input$blast_table_rows_selected))
    y = z[ x ]
    points(x, y, col='navyblue', pch=1, cex=2.25)
  }

  pos <- c(rep(3, length(gp))[-length(gp)],2)
  text(gp, z[gp],
       labels=paste0("Nhit=", gp, ", cutoff=", round(z[gp])),
       col="black", pos=3, cex=1)

  legend("bottomleft", c("Selected hits", "Above cutoff", "Below cutoff", "Cutoff"),
         pt.bg = c("blue", "red", "grey50", "red"), col=c("blue", rep("grey10", 2), "red"), bg="grey90", lty=c(0, 0, 0, 2), pch=c(1,21,21,NA))


})

## Currently not used
output$blast_plot2 <- renderChart2({
  ##blast <- run_blast()
  blast <- rv$blast
  cut <- set_cutoff(blast, rv$cutoff)
  blast$mlog.evalue <- blast$score

  gp <- cut$gp.inds
  cutoff <- cut$cutoff
  grps <- cut$grps
  z <- blast$score

  ## generate a dataframe: pc$z + pdbids + group
  data <- data.frame(
    Bitscore = z,
    id = blast$acc,
    group = sapply(grps, function(x) if(x==1) "Above" else "Below cutoff"),
    Hits = 1:length(z),
    desc = blast$desc
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
      return '<b>pdb:</b> ' + e.point.id + '</br>' + e.point.desc
    } !#")
    p1$setTemplate(afterScript = '<script>
      var css = document.createElement("style");
      css.type = "text/css";
      css.innerHTML = ".nv-x .nv-axislabel { font-size: 20px; }";
      document.body.appendChild(css);
      css = document.createElement("style");
      css.innerHTML = ".nv-y .nv-axislabel { font-size: 20px; }";
      document.body.appendChild(css);
      document.getElementById("blast_plot2").focus();
    </script>'
                   )

  return(p1)
})


get_blasttable <- reactive({
  message("fetching blast table")

  progress <- shiny::Progress$new(session, min=1, max=5)
  progress$set(message = 'Annotating results',
               detail = 'Please wait ...')
  progress$set(value = 1)
  
  ptm <- proc.time()

  if(input$input_type != "multipdb") {

    blast <- rv$blast
    grps <- set_cutoff(blast, cutoff=rv$cutoff)$grps

    acc <- blast$acc
    anno <- get_annotation(acc)
    progress$set(value = 2)

    ## clean up if length does not correspond
    if(length(acc) != length(anno$acc)) {
      print(paste("blast: ", length(acc)))
      print(paste("anno: ", length(anno$acc)))

      if(length(acc) > length(anno$acc)) {
        inds <- acc %in% anno$acc
        missed <- blast$acc[!inds]
        blast <- blast[inds,, drop=FALSE]

        print("could not find annotation for pdb ids")
        print(missed)

        print(paste("corrected blast: ", length(blast$acc)))
        print(paste("anno: ", length(anno$acc)))
        acc <- blast$acc
      }
    }

    anno$score <- blast$score
  }
  else {
    acc <- unique(trim(unlist(strsplit(input$pdb_codes, ","))))
    acc <- acc[acc!=""]

    if(!length(acc) > 0)
      return()

    acc <- format_pdbids(acc)
    anno <- get_annotation(acc, use_chain=FALSE)

    inds <- unlist(sapply(acc, grep, anno$acc))
    anno <- anno[inds, ]
    acc <- anno$acc

    anno$score <- rep(0, length(acc))

    hits <- NULL
    grps <- NULL
    progress$set(value = 2)
  }

  anno$struct_nr <- 1:nrow(anno)

  if(!is.null(grps)) {
    col <- sapply(grps[1:nrow(anno)], function(x) { if(x==1) 'red' else 'black' })

    #if(length(input$selected_pdbids)>0) {
    #  inds <- unlist(lapply(input$selected_pdbids, grep, anno$acc))
    #  col[inds] <- "green"
    #}

    #if(!is.null(rv$limit_hits)) {
    #  limit <- as.numeric(rv$limit_hits)
    #  if(limit > 0)
    #    col[1:limit] <- "green"
    #}
    anno$acc <- paste0(
      "<span class=\"id_", col, "\" style=\"color:",
      col,
      "; font-size:large\">&#x25CF;</span>",
      "&nbsp;<a href=\"", "http://pdb.org/pdb/explore/explore.do?structureId=",
      substr(anno$acc, 1, 4), "\" target=\"_blank\">", anno$acc, "</a>"
      )
  }
  else {
#    anno$acc <- paste0("<a href=\"", "http://pdb.org/pdb/explore/explore.do?structureId=",
#                         substr(anno$acc, 1, 4), "\" target=\"_blank\">", anno$acc, "</a>")
    anno$acc <- sapply(anno$acc, function(x) { as.character(tags$a(href = paste0('http://pdb.org/pdb/explore/explore.do?structureId=', substr(x, 1, 4)), target = '_blank', x)) })
  }

  progress$set(value = 3)

  ## http://xray.bmc.uu.se/hicup/ATP/
  anno$ligandId <- unlist(lapply(anno$ligandId,
         function(x) {
           if(is.na(x))
             return("")

           spl <- unlist(strsplit(x, ","))
           if(length(spl) > 0) {
             return(paste(
                         sapply(spl, function(x) { as.character(tags$a(href = paste0("http://www.rcsb.org/pdb/ligand/ligandsummary.do?hetId=", x), target = '_blank', x))}),
                          collapse=", "))
           }
           else {
             return("")
           }

         }))

  progress$set(value = 4)
  rownames(anno) <- NULL

  show.cols <- c("acc", "compound", "source", "ligandId", "score", "experimentalTechnique", "resolution", "ligandName")

  col.inds <- sapply(show.cols, grep, colnames(anno))

  message("time used for fetching blast table:")
  print(proc.time() - ptm)
  progress$set(value = 5)

  gc()
  progress$close()
  
  # re-ordering changed below to ID, Score .. for datatable
  return(anno[, col.inds[c(1,5,2,3,4, 6:8)]])
})


  output$blast_table <- renderDataTable({
      limit <- as.numeric(input$limit_hits)
      DT::datatable(get_blasttable(), extensions = c('Scroller', 'ColVis'),## 'TableTools'),
              class = 'compact stripe cell-border',   ## Remove 'compact' to add padding
              escape = FALSE,
              ##- Note. Should have 'PFAM' and 'Authors' here also...
              colnames = c("ID", "BitScore", "Name (PDB Title)", "Species",
                           "Ligands", "Method", "Resolution (A)", "Ligand Names"),

              selection =
                  if(input$input_type != 'multipdb') {
                      list(mode = 'multiple',
                          selected =
                              if(length(limit) == 0) {
                                  as.character(1:5)
                              } else {
                                  as.character(1:limit)
                              }
                      )
                  } else {
                      "none"
                  },

              options =
                list(dom='C<"clear">frtiS', ##-- with ColVis 'C' and Scroler 'S' options
                colVis = list(exclude = c(0, 1), activate = 'click', buttonText="Show/Hide Columns"),

                ##dom = "frtiS",               ##--- ORIG
                ##   dom = 'T<"clear">lfrtip', ##-- with tableTools option
                ##   tableTools = list(sSwfPath = copySWF()),

                scrollY = 400,                 ##---- Scroler options (height)
                scrollX = '100%',
                scrollXInner = '100%',
                scrollCollapse = TRUE,
                deferRender = TRUE,

                autoWidth = FALSE,

                columnDefs = list(                          ##-  Set cols behavior
                  list( width = '5%', targets = c(0) ),     ## nums
                  list( width = '10%', targets = c(1, 2) ), ## ID and BitScore
                  list( width = '40%', targets = c(3) ),    ## Name
                  list( width = '20%', targets = c(4) ),    ## Species
                  list( width = '15%', targets = c(5) ),    ## Ligands

                  list( visible=FALSE, targets = c(6,7,8) ), ## Others.. R-free, Authors, PubMed...
                  list( orderable = 'false', targets = c(0,8) )
                ),

                initComplete = JS(
                'function(settings) {',
                'console.log("search datatable loaded");',
                '}'
                )
              )
      )
  }) ## End renderDT


output$cutoff_slider <- renderUI({
  ##blast <- run_blast()
  blast <- rv$blast
  cutoff <- rv$cutoff

  sliderInput("cutoff", "Adjust inclusion BitScore cutoff:",
              min = floor(min(blast$score)), max = floor(max(blast$score)), value = cutoff)
})

output$hits_slider <- renderUI({
  hits <- filter_hits()

  sliderInput("limit_hits", "Limit total number of included structures:",
              min = 1, max = length(hits$hits), value = rv$limit_hits, step=1)
})

