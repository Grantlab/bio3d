library(staticdocs)
build_package("bio3d", base_path="html", examples=TRUE)

save.image("tmp_build_package.RData")
q("no")

# run explicitly
package_path="../bio3d"
base_path="html"
#examples=FALSE
examples=TRUE
devtools:::load_all(package_path)
package <- staticdocs:::package_info(package_path, base_path, examples)
if (!file.exists(package$base_path)) 
    dir.create(package$base_path)
staticdocs:::copy_bootstrap(base_path)
#package$topics <- staticdocs:::build_topics(package)
    index <- package$rd_index
    paths <- file.path(package$base_path, index$file_out)
    index$title <- ""
    index$in_index <- TRUE
    for (i in seq_along(index$name)) {
        message("Generating ", basename(paths[[i]]))
        rd <- package$rd[[i]]
        html <- to_html(rd, env = new.env(parent = globalenv()), 
            topic = str_replace(basename(paths[[i]]), "\\.html$", 
                ""), package = package)
        html$pagetitle <- html$name
        html$package <- package[c("package", "version")]
        render_page(package, "topic", html, paths[[i]])
        graphics.off()
        if ("internal" %in% html$keywords) {
            index$in_index[i] <- FALSE
        }
        index$title[i] <- html$title
    }
    index

package$vignettes <- staticdocs:::build_vignettes(package)
package$demos <- staticdocs:::build_demos(package)
package$readme <- staticdocs:::readme(package)
staticdocs:::build_index(package)
if (interactive()) {
    browseURL(normalizePath(file.path(base_path, "index.html")))
}
invisible(TRUE)
save.image("tmp_build_package.RData")
