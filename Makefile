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
