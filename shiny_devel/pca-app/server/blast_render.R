output$pdbWebGL  <- renderWebGL({
  pdb <- get_pdb()
  view.pdb(pdb, as=input$view_inpdb_as, col=input$view_inpdb_col)
})

output$blast_plot <- renderUI({
  blast <- run_blast()

  if(length(blast$acc) > 100) {
    plotOutput("blast_plot1")
  }
  else {
    showOutput("blast_plot2", "nvd3")

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
  }

})

### Blast plot on front page
output$blast_plot1 <- renderPlot({
  blast <- run_blast()
  cut <- set_cutoff(blast, rv$cutoff)

  gp <- cut$gp.inds
  cutoff <- cut$cutoff
  grps <- cut$grps
  z <- blast$score

  #col <- sapply(grps, function(x) if(x==1) 'red' else if(x==2) 'grey50')
  #if(length(input$selected_pdbids)>0) {
  #  inds <- unlist(lapply(input$selected_pdbids, grep, blast$acc))
  #  col[inds] <- "green"
  #}
  #else {
  #  limit <- as.numeric(rv$limit_hits)
  #  if(limit > 0)
  #    col[1:limit] <- "green"
  #}

  ##if(!length(input$blast_table_rows_selected)>0) {
  #if(!length(input$hits_in)>0) {
  #  limit <- as.numeric(rv$limit_hits)
  #  if(limit > 0)
  #    col[1:limit] <- "green"
  #}
  #else {
  #  col[ as.numeric(input$hits_in) ] <- "green"
  #}

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
  pch = rep(19, length(col.bg))
  pch[col.bg=='grey50'] <- 21
  col.pch = rep('red3', length(pch))
  col.pch[pch==21] <- 'black'

  par(mar=c(6, 4, 0, 0))
  plot(z, xlab="", ylab="Bitscore of Alignment to Input", bg=col.bg, col=col.pch, pch=pch, cex=1.1)
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
         pt.bg = c("green", "red", "grey50", "red"), col=c(rep("grey10", 3), "red"), bg="grey90", lty=c(0, 0, 0, 2), pch=c(21,21,21,NA))

})


output$blast_plot2 <- renderChart2({
  blast <- run_blast()
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

  if(input$input_type != "multipdb") {
    blast <- run_blast()
    grps <- set_cutoff(blast, cutoff=rv$cutoff)$grps

    acc <- blast$acc
    anno <- get_annotation(acc)

    ## clean up if length does not correspond
    if(length(acc) != length(anno$acc)) {
      print(paste("blast: ", length(acc)))
      print(paste("anno: ", length(anno$acc)))

      if(length(acc) > length(anno$acc)) {
        inds <- acc %in% anno$acc
        blast <- blast[inds,, drop=FALSE]
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
  }

  anno$struct_nr <- 1:nrow(anno)

  if(!is.null(grps)) {
    col <- sapply(grps[1:nrow(anno)], function(x) { if(x==1) 'red' else 'black' })

    if(length(input$selected_pdbids)>0) {
      inds <- unlist(lapply(input$selected_pdbids, grep, anno$acc))
      col[inds] <- "green"
    }

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
    anno$acc <- paste0("<a href=\"", "http://pdb.org/pdb/explore/explore.do?structureId=",
                         substr(anno$acc, 1, 4), "\" target=\"_blank\">", anno$acc, "</a>")
  }

  rownames(anno) <- NULL
  show.cols <- c("acc", "compound", "source", "ligandId", "score")
  col.inds <- sapply(show.cols, grep, colnames(anno))

  # re-ordering changed below to ID, Score .. for datatable
  return(anno[, col.inds[c(1,5,2,3,4)]])
})


  output$blast_table <- renderDataTable({
      x <- get_blasttable()
      limit <- as.numeric(input$limit_hits)
      DT::datatable(x, extensions = 'Scroller', escape = FALSE,
              colnames = c("ID", "Score", "Name", "Species", "Ligands"),
              selection = list(mode = 'multiple',
                               selected = if(length(limit) == 0) as.character(1:5) else as.character(1:limit)),
              options = list(
                deferRender = TRUE,
                dom = "frtiS",
                scrollY = 400,
                scrollCollapse = TRUE,
                autoWidth = FALSE,
                columnDefs = list(
                  list( width = '5%', targets = c(0) ),
                  list( width = '10%', targets = c(1, 2) ),
                  list( width = '40%', targets = c(3) ),
                  list( width = '20%', targets = c(4) ),
                  list( width = '15%', targets = c(5) )
                ),
                initComplete = JS(
                'function(settings) {',
                'document.getElementById("blast_plot").scrollIntoView();',
                'console.log("initComplete complete");',
                '}'
                )
             )
      )
  })

## checkbox
output$pdb_chains <- renderUI({
  chains <- get_chainids()

  selectInput("chainId", "Limit to chain ID:",
              choices = chains, selected = chains[1], multiple=FALSE)
})

output$cutoff_slider <- renderUI({
  blast <- run_blast()
  cutoff <- rv$cutoff

  sliderInput("cutoff", "Adjust inclusion bitscore cutoff:",
              min = floor(min(blast$score)), max = floor(max(blast$score)), value = cutoff)
})

output$hits_slider <- renderUI({
  hits <- filter_hits()

  sliderInput("limit_hits", "Limit total number of included structures:",
              min = 1, max = length(hits$hits), value = rv$limit_hits, step=1)
})

