library(staticdocs)
build_package("~/bio3d/ver_devel/bio3d", base_path="~/bio3d/ver_devel/util/html", examples=TRUE)

save.image("tmp_build_package.RData")
q("no")

# run explicitly
package="~/bio3d/ver_devel/bio3d"
base_path="~/bio3d/ver_devel/util/html"
#examples=FALSE
examples=TRUE
devtools:::load_all(package)
package <- staticdocs:::package_info(package, base_path, examples)
if (!file.exists(package$base_path)) 
    dir.create(package$base_path)
staticdocs:::copy_bootstrap(base_path)
package$topics <- staticdocs:::build_topics(package)
package$vignettes <- staticdocs:::build_vignettes(package)
package$demos <- staticdocs:::build_demos(package)
package$readme <- staticdocs:::readme(package)
staticdocs:::build_index(package)
if (interactive()) {
    browseURL(normalizePath(file.path(base_path, "index.html")))
}
invisible(TRUE)
save.image("tmp_build_package.RData")
