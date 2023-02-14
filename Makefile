SHELL=/bin/bash
.PHONY: all docker
all:
	./gradlew jvarkit

docker: Dockerfile
	# build
	docker build -t `grep hard $<  | tr " " "\n" | grep hard -A 1 | tail -1 | awk '{printf("lindenb/jvarkit:%s",$$1);}'` .
	# list images
	docker images


docker-clean:
	docker image rm `docker images  |grep 'lindenb/jvarkit' | tr -s " "  | cut -d ' ' -f 3`
