HOMEPATH=$(PWD)
MODLOC=/usr2/postdoc/piskulic/privatemodules

pair_dist:
	@echo "Making pair distribution code"
	cp modulefiles/pair-distribution.lua $(MODLOC)
	@echo "prepend_path('PATH', \"$(HOMEPATH)/bin\")" >> $(MODLOC)/pair-distribution.lua
	mkdir -p bin/
	touch bin/test
	rm bin/*
	ln -s $(HOMEPATH)/lipid-rdf-code/setup_rdf.py $(HOMEPATH)/bin/
	ln -s $(HOMEPATH)/lipid-rdf-code/rdf_main.py $(HOMEPATH)/bin/
	ln -s $(HOMEPATH)/lipid-rdf-code/dump_split.py $(HOMEPATH)/bin/
	chmod 777 bin/*
