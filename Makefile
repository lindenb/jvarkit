SHELL=/bin/bash
.PHONY: all docker
all:
	./gradlew jvarkit

docker: Dockerfile
	# build
	docker build -t `grep -m1 "JVARKIT_VERSION=" $<  | cut -d '=' -f2 | awk '{printf("lindenb/jvarkit:%s",$$1);}'` .
	# list images
	docker images


docker-clean:
	docker image rm `docker images  |grep 'lindenb/jvarkit' | tr -s " "  | cut -d ' ' -f 3`


schemas : dtd-pubmed dtd-blast dtd-insdseq dtd-gb xsd-uniprot

dtd-pubmed:
	./xjc -nv -d src/main/java -p gov.nih.nlm.ncbi.pubmed -dtd  "https://www.ncbi.nlm.nih.gov/corehtml/query/DTD/pubmed_100101.dtd"

dtd-blast:
	./xjc -nv -d src/main/java -p gov.nih.nlm.ncbi.blast -dtd  "https://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd"

dtd-insdseq:
	./xjc -nv -d src/main/java -p gov.nih.nlm.ncbi.insdseq -dtd  "https://www.ncbi.nlm.nih.gov/dtd/INSD_INSDSeq.dtd"

dtd-gb:
	./xjc -nv -d src/main/java -p gov.nih.nlm.ncbi.gb -dtd   "https://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.dtd"

xsd-uniprot:
	./xjc -nv -d src/main/java -p org.uniprot "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot.xsd"
