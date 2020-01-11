outfile:
	@echo $(basename $(INFILE)).png

build:
	plantuml $(INFILE)

view:
	feh -R 2 $(OUTFILE) 2>/dev/null

refresh:
