

rawlines = readLines("tmp.md")
newlines <- rawlines

fig.start <- grep("!\\[", rawlines)
fig.end <- grep(".png)", rawlines)


for ( i in 1:length(fig.start)) {
  txt <- paste(rawlines[fig.start[i]:fig.end[i]], collapse=" ")

  sep <- unlist(strsplit(txt, "\\]\\("))
  caption <- paste("<strong>Figure ", i, "</strong>: ",
                   substr(sep[1], 3, nchar(sep[1])), sep="")
  
  src <- substr(sep[2], 1, nchar(sep[2])-1)
  
  html <- paste('<figure>', 
                '<img src="', src, '">',
                '<figcaption>', caption, '<figcaptionend>', '<figureend>', sep="")

  newlines[fig.start[i]:fig.end[i]] <- ""
  newlines[fig.start[i]] <- html
  
}

writeLines(newlines, "new.md")
