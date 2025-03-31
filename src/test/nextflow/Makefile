NEXTFLOW0=$(realpath ${HOME}/src/nextflow/nextflow)
ifeq ($(NEXTFLOW0),'')
NEXTFLOW=nextflow
else
NEXTFLOW=$(NEXTFLOW0)
endif
all: main.nf nextflow.config ../../../dist/jvarkit.jar
	$(NEXTFLOW) run -work-dir work -resume $< --testDir "$(realpath ../../test/resources)" --jvarkit "$(realpath ../../../dist/jvarkit.jar)"
	
clean:
	rm -rf work .nextflow.log .nextflow*
