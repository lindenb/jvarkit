
SELF=$(realpath $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST)))
EXEC=java -Djava.io.tmpdir=. -jar ../dist/jvarkit.jar
RSRC=../src/test/resources
all: $(shell grep '^test' "$(SELF)" | cut -d ':' -f1)
	@echo "COMPLETED $(words $^) TESTS"

test.dict2xml:
	$(EXEC) dict2xml $(RSRC)/S1.bam
	$(EXEC) dict2xml $(RSRC)/S1.vcf.gz
	echo $(RSRC)/S1.bam  | $(EXEC) dict2xml
	$(EXEC) dict2xml "https://raw.githubusercontent.com/VertebrateResequencing/vr-codebase/653b76ad0590e2ccaf9dd9fd6ad4596072045c02/t/data/bam_exome_qc_bait.interval_list"

test.vcf2intervals:
	$(EXEC) vcf2intervals --bed -D 100 $(RSRC)/rotavirus_rf.vcf.gz
	$(EXEC) vcf2intervals --bed -N 10 $(RSRC)/rotavirus_rf.vcf.gz
	$(EXEC) vcf2intervals -D 100 $(RSRC)/rotavirus_rf.vcf.gz > jeter.interval_list && \
		$(EXEC) dict2xml jeter.interval_list && \
		rm jeter.interval_list
	$(EXEC) vcf2intervals -N 10 $(RSRC)/rotavirus_rf.vcf.gz > jeter.interval_list && \
		$(EXEC) dict2xml jeter.interval_list && \
		rm jeter.interval_list

test.vcffilterjdk:
	$(EXEC) vcffilterjdk -e 'return variant.getStart()%2==0;' $(RSRC)/rotavirus_rf.vcf.gz

test.vcfpolyx:
	$(EXEC) vcfpolyx -R $(RSRC)/rotavirus_rf.fa $(RSRC)/rotavirus_rf.vcf.gz
	$(EXEC) vcfpolyx -R $(RSRC)/rotavirus_rf.fa -n 2 --tag POLYX $(RSRC)/rotavirus_rf.vcf.gz

test.vcftail:
	$(EXEC) vcftail -n 2 $(RSRC)/rotavirus_rf.vcf.gz | grep -v '^#' | wc -l | grep -wF 2

test.vcfhead:
	$(EXEC) vcfhead -n 2 $(RSRC)/rotavirus_rf.vcf.gz | grep -v '^#' | wc -l | grep -wF 2

test.vcfrebase:
	$(EXEC) vcfrebase -R $(RSRC)/rotavirus_rf.fa $(RSRC)/rotavirus_rf.vcf.gz

test.wescnvsvg:
	$(EXEC) wescnvsvg -R $(RSRC)/rotavirus_rf.fa  -r 'RF01:1-1000'  $(RSRC)/S1.bam $(RSRC)/S2.bam $(RSRC)/S3.bam

test.wgscoverageplotter:
	$(EXEC) wgscoverageplotter -R $(RSRC)/rotavirus_rf.fa $(RSRC)/S1.bam
