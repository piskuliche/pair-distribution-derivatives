HOMEPATH=$(PWD)
MODLOC="/usr2/postdoc/piskulic/privatemodules"

pair_dist:
	@echo "Making pair distribution code"
	cp src/modulefiles/pair-distribution.lua $(MODLOC)
	@echo "prepend_path('PATH',$(HOMEPATH)/bin)" >> $(MODLOC)/pair-distribution.lua
	mkdir -p bin/
	touch bin/test
	rm bin/*
	ln -s $(HOMEPATH)/src/setup_rdf.py $(HOMEPATH)/bin/
	ln -s $(HOMEPATH)/src/rdf_main.py $(HOMEPATH)/bin/
	chmod 777 bin/*