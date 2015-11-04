.PHONY:all tests
REF=${htsjdk.home}/testdata/htsjdk/samtools/reference/Homo_sapiens_assembly18.trimmed.fasta
OUTPUT=tests_results.txt
define compile
	(cd .. && make $(1)) || (echo "$(1)	Compilation Failed" >> ${OUTPUT} && exit 999)
endef


TEST_APPS = vcffilterjs bamtile addlinearindextobed bamcmpcoverage allelefreqcalc vcfhead vcftail

all:tests

tests: ${TEST_APPS}

vcffilterjs:
	$(call compile,$@) && java -jar ${dist.dir}/$@.jar -e 'variant.getStart()%2==0' ${htsjdk.home}/testdata/htsjdk/variant/ILLUMINA.wex.broad_phase2_baseline.20111114.both.exome.genotypes.1000.vcf > /dev/null && echo "OK: $@" > $(addsuffix .ok,$@)

allelefreqcalc: 
	$(call compile,$@) && java -jar ${dist.dir}/$@.jar ${htsjdk.home}/testdata/htsjdk/variant/ILLUMINA.wex.broad_phase2_baseline.20111114.both.exome.genotypes.1000.vcf > /dev/null && echo "OK: $@" > $(addsuffix .ok,$@)

addlinearindextobed bamcmpcoverage:	
	$(call compile,$@) && java -jar ${dist.dir}/$@.jar -h > /dev/null && echo "OK: $@" > $(addsuffix .ok,$@)

bamtile:
	$(call compile,$@) && java -jar ${dist.dir}/$@.jar ${htsjdk.home}/testdata/htsjdk/samtools/compressed.bam  && echo "OK: $@" > $(addsuffix .ok,$@)

vcfhead:
	$(call compile,$@) && java -jar ${dist.dir}/$@.jar -n 2 ${htsjdk.home}/testdata/htsjdk/variant/ILLUMINA.wex.broad_phase2_baseline.20111114.both.exome.genotypes.1000.vcf > /dev/null && echo "OK: $@" > $(addsuffix .ok,$@)

vcftail:
	$(call compile,$@) && java -jar ${dist.dir}/$@.jar -n 2 ${htsjdk.home}/testdata/htsjdk/variant/ILLUMINA.wex.broad_phase2_baseline.20111114.both.exome.genotypes.1000.vcf > /dev/null && echo "OK: $@" > $(addsuffix .ok,$@)




vcftreepack:
	$(call compile,$@) && echo "$@"

vcftrio:
	$(call compile,$@) && echo "$@"

vcfucsc:
	$(call compile,$@) && echo "$@"

vcfvcf:
	$(call compile,$@) && echo "$@"

vcfviewgui:
	$(call compile,$@) && echo "$@"

vcfviewgui.exe:
	$(call compile,$@) && echo "$@"

worldmapgenome:
	$(call compile,$@) && echo "$@"

