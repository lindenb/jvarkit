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


schemas : dtd-pubmed dtd.blast

dtd-pubmed:
	./xjc -nv -d src/main/java -p com.github.lindenb.jvarkit.ncbi.schema.pubmed -dtd  "https://www.ncbi.nlm.nih.gov/corehtml/query/DTD/pubmed_100101.dtd"

dtd-blast:
	./xjc -nv -d src/main/java -p com.github.lindenb.jvarkit.ncbi.schema.blast -dtd  "https://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd"
