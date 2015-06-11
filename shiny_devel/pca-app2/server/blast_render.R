output$pdbWebGL  <- renderWebGL({
  pdb <- get_pdb()
  view.pdb(pdb, as="overview", col="sse")
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
  print(paste0("dimZ: ",dim(z)))

  col.bg <- sapply(grps, function(x) if(x==1) 'red3' else if(x==2) 'grey50')

  pch = rep(19, length(col.bg))
  pch[col.bg=='grey50'] <- 21
  col.pch = rep('red3', length(pch))
  col.pch[pch==21] <- 'black'
  plot(z, xlab="", ylab="Bitscore", bg=col.bg, col=col.pch, pch=pch, cex=1.1)
  abline(v=gp, col="gray50", lty=3)
  abline(h=cutoff, col="red", lty=2)
  print(as.numeric(input$blast_table_rows_selected))
  if(!length(input$blast_table_rows_selected)>0) {
    limit <- sort(as.numeric(input$limit_hits))
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
         pt.bg = c(NA, "red3", "grey50", "red"), col=c('navyblue', 'red3', 'grey10', 'red'), bg="grey90", lty=c(0, 0, 0, 2), pch=c(1,19,21,NA))

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


get_blasttable <- reactive({
  if(input$input_type != "multipdb") {
    blast <- run_blast()
    hits <- filter_hits()
    grps <- set_cutoff(blast, cutoff=input$cutoff)$grps

    ## provide additional 5 non-checked table rows
    checked.inds <- which(hits$hits)
    unchecked.inds <- which(!hits$hits)

    n <- length(blast$acc)
    m <- n - sum(hits$hits)
    #if (m > 5)
    #  m <- 5

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

  anno$score <- hits$score
  anno$struct_nr <- 1:nrow(anno)

  if(!is.null(grps)) {
    col <- sapply(grps[1:nrow(anno)], function(x) { if(x==1) 'red' else 'black' })

    #if(!is.null(input$limit_hits)) {
    #  limit <- as.numeric(input$limit_hits)
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

  checkbox <- paste0("<input type=\"checkbox\" name=\"pdb_ids\" value=\"", hits$acc, "\" ",  checked, ">")
  anno$check <- checkbox

  rownames(anno) <- NULL
  show.cols <- c("acc", "compound", "source", "ligandId", "score")
  col.inds <- sapply(show.cols, grep, colnames(anno))

  #print(col.inds)
  #print(head(anno[, col.inds]))
  return(anno[, col.inds])
})

output$blast_table <- renderDataTable({
    DT::datatable(get_blasttable(), extensions = 'Scroller', escape = FALSE,
            colnames = c("ID", "Name", "Species", "Ligands", "Score"),
            options = list(
              deferRender = TRUE,
              dom = "frtiS",
              scrollY = 400,
              scrollCollapse = TRUE,
              initComplete = JS(
                'function(settings, json) {',
                '//var table = $("#DataTables_Table_1").dataTable();',
                'var nodes = this.fnGetNodes(); //fnGetDisplayNodes',
                paste0('for(var i =0; i < ', if(length(input$limit_hits)==0) 5 else as.numeric(input$limit_hits), '; i++) {'),
                '$(nodes[i]).addClass("selected");',
                '}',
                '$(nodes[0]).click();',
                '$(nodes[0]).click();',
                '}'
              )
           )
    )
})




## checkbox
output$pdb_chains <- renderUI({
  chains <- get_chainids()
  ##radioButtons("chainId", label="Choose chain ID:",
  ##             choices=chains, inline=TRUE)

  selectInput("chainId", "Limit to chain ID:",
              choices = chains, selected = chains[1], multiple=FALSE)
  
})


output$resetable_cutoff_slider <- renderUI({
  reset <- input$reset_cutoff

  blast <- run_blast()
  cutoff <- set_cutoff(blast, rv$cutoff)$cutoff

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
  hits <- set_cutoff(blast, cutoff=rv$cutoff)

  sliderInput("limit_hits", "Limit hits:",
              min = 1, max = length(hits$inds), value = 5, step=1)
})

