ROOTDIR= .
include $(ROOTDIR)/Makefile.base



SUBDIRS = utility matrices MLGkernel MLGkernel/swig # tests

.PHONY: all objects tests clean $(SUBDIRS)

intro:
	@echo; echo "\033[1mCompiling MLGkernel system... \033[00m"; echo 
ifdef EIGENDIR
	@echo "Linear algebra with Eigen:            enabled"
else
	@echo "Linear algebra with Eigen:           disabled"
endif


objects:
	@for dir in $(SUBDIRS); do echo; echo "\033[1m *** Making objects in $$dir directory *** \033[00m"; $(MAKE) -C $$dir objects; done

all: intro
	@for dir in $(SUBDIRS); do echo; echo "\033[1m *** Making all in $$dir directory *** \033[00m"; $(MAKE) -C $$dir all; done; echo "\n\n" 

clean:
	@for dir in $(SUBDIRS); do $(MAKE) -C $$dir clean; done


# declare dependencies among subdirectories
filetypes: utility
matrices: filetypes

anew: clean all 
