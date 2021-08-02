#!usr/bin/make -f
# All commands are run as R functions rather than shell commands so that it will work easily on any Windows machine, even if the Windows machine isn't properly set up with all the right tools

all: README.md

clean:
	Rscript -e 'suppressWarnings(file.remove("README.md", "vignettes/vignette.Rmd"))'

.PHONY: all clean
.DELETE_ON_ERROR:
.SECONDARY:

README.md : vignettes/vignette.Rmd
#	echo "Rendering the overview vignette"
	Rscript -e 'rmarkdown::render("vignettes/vignette.Rmd", output_format = "md_document")'
#	echo "Correcting image paths"
#	sed -i -- 's,../inst,inst,g' vignettes/vignette.md
	Rscript -e 'file <- gsub("../inst", "inst", readLines("vignettes/vignette.md")); writeLines(file, "vignettes/vignette.md")'
#	echo "Copying output to README.md"
#	cp vignettes/vignette.md README.md
	Rscript -e 'file.copy("vignettes/vignette.md", "README.md", overwrite = TRUE)'
	Rscript -e 'suppressWarnings(file.remove("vignettes/vignette.md"))'
