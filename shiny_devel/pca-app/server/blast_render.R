output$pdbWebGL  <- renderWebGL({
  pdb <- get_pdb()
  view.pdb(pdb)
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
  hits <- set_cutoff(blast, input$cutoff)
  cutoff <- hits$cutoff
  gp <- hits$gp.inds
  grps <- hits$grps
  z <- blast$score

  col <- sapply(grps, function(x) if(x==1) 'red' else if(x==2) 'black')
  l <- as.numeric(input$limit_hits)
  col[1:l] <- "blue"

  plot(z, xlab="", ylab="Bitscore", bg=col, col="grey10", pch=21, cex=1.1)
  abline(v=gp, col="gray50", lty=3)
  abline(h=cutoff, col="gray50", lty=3)

  pos <- c(rep(3, length(gp))[-length(gp)],2)
  text(gp, z[gp],
       labels=paste0("Nhit=", gp, ", cutoff=", round(z[gp])),
       col="black", pos=3, cex=1)

  legend("bottomleft", c("Selected hits", "Above cutoff", "Below cutoff"),
         col=c("blue", "red", "black"), pch=16)

})


output$blast_plot2 <- renderChart2({
  blast <- run_blast()
  hits <- set_cutoff(blast, input$cutoff)
  blast$mlog.evalue <- blast$score

  cutoff <- hits$cutoff
  gp <- hits$gp.hits
  gps <- hits$grps
  z <- blast$score

  ## generate a dataframe: pc$z + pdbids + group
  data <- data.frame(
    Bitscore = z,
    id = blast$acc,
    group = sapply(gps, function(x) if(x==1) "Above" else "Below cutoff"),
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


### Table of annotated BLAST results
output$blast_table <- renderDataTable({
  message("rendering blast data table")
  if(input$input_type != "multipdb") {
    blast <- run_blast()
    hits <- filter_hits()
    grps <- set_cutoff(blast, cutoff=input$cutoff)$grps

    ## provide additional 5 non-checked table rows
    checked.inds <- which(hits$hits)
    unchecked.inds <- which(!hits$hits)

    n <- length(blast$acc)
    m <- n - sum(hits$hits)
    if (m > 5)
      m <- 5

    unchecked.inds <- unchecked.inds[1:m]
    show.inds <- c(checked.inds, unchecked.inds)

    hits <- blast[show.inds,, drop=FALSE]
    acc <- hits$acc
    anno <- get_annotation(acc)

    checked <- rep("CHECKED", length(acc))
    checked[unchecked.inds] <- ""
  }
  else {
    acc <- toupper(unique(trim(unlist(strsplit(input$pdb_codes, ",")))))
    acc <- acc[acc!=""]
    anno <- get_annotation(acc, use_chain=FALSE)
    inds <- unlist(sapply(acc, grep, anno$acc))
    anno <- anno[inds, ]
    acc <- anno$acc

    hits <- NULL
    grps <- NULL
    hits$acc <- acc
    
    hits$score <- rep(0, length(acc))
    checked <- rep("CHECKED", length(acc))
  }

  if(!is.null(grps) & !is.null(grps)) {
    col <- sapply(grps[1:nrow(anno)], function(x) { if(x==1) 'red' else 'black' })
    l <- as.numeric(input$limit_hits)
    col[1:l] <- "blue"
    
    anno$pdbId <- paste0(
      "<span style=\"color:",
      col,
      "; font-size:large\">&#x25CF;</span>",
      "&nbsp;<a href=\"", "http://pdb.org/pdb/explore/explore.do?structureId=",
      substr(anno$acc, 1, 4), "\" target=\"_blank\">", anno$acc, "</a>"
    )
  } else {
      anno$pdbId <- paste0("<a href=\"", "http://pdb.org/pdb/explore/explore.do?structureId=",
        substr(anno$acc, 1, 4), "\" target=\"_blank\">", anno$acc, "</a>")
  }
  anno$score <- hits$score
  anno$id <- 1:nrow(anno)

  checkbox <- paste0("<input type=\"checkbox\" name=\"pdb_ids\" value=\"", hits$acc, "\"",  checked, ">")
  anno$check <- checkbox

  return(anno[, c("id", "check", "pdbId", "compound", "source", "ligandId", "chainLength", "score")])
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
  cutoff <- set_cutoff(blast, input$cutoff)$cutoff

  sliderInput("cutoff", "Adjust cutoff:",
              min = floor(min(blast$score)), max = floor(max(blast$score)), value = cutoff)
})

output$resetable_cutoff_slider <- renderUI({
  reset <- input$reset_cutoff

  blast <- run_blast()
  hits <- set_cutoff(blast, cutoff=NULL)
  cutoff <- hits$cutoff

  sliderInput("cutoff", "Adjust cutoff:",
              min = floor(min(blast$score)), max = floor(max(blast$score)), value = cutoff)
})

output$hits_slider <- renderUI({
  blast <- run_blast()
  hits <- set_cutoff(blast, cutoff=input$cutoff)

  sliderInput("limit_hits", "Limit hits:",
              min = 1, max = length(hits$inds), value = 5, step=1)
})

