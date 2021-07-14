FROM alpine:latest AS build-base
RUN apk add build-base util-linux linux-headers g++ bash perl make cmake git bison flex openjdk11 apache-ant
ENV JAVA_HOME /usr/lib/jvm/java-11-openjdk
ENV PATH="${PATH}:${JAVA_HOME}/bin"

FROM build-base AS build
ARG VDB_BRANCH=engineering
ARG NGS_BRANCH=engineering
ARG NGT_BRANCH=tax
RUN git clone --branch ${NGS_BRANCH} --depth 1 https://github.com/ncbi/ngs.git
RUN git clone --branch ${VDB_BRANCH} --depth 1 https://github.com/ncbi/ncbi-vdb.git
RUN git clone --branch ${NGT_BRANCH} https://github.com/ncbi/ngs-tools.git
ARG BUILD_STYLE=--without-debug
WORKDIR /ngs
RUN ./configure ${BUILD_STYLE} && make -s -C ngs-sdk
WORKDIR /ncbi-vdb
RUN sed -e '/gnu\/libc-version.h/ s/^/\/\//' -e '/gnu_get_libc_version/s/^/\/\//' -i libs/kns/manager.c
RUN sed -e '/\#if _DEBUGGING && LINUX/ s/$/ \&\& 0/' -i interfaces/kfc/defs.h
RUN ./configure ${BUILD_STYLE} && make -s
WORKDIR /ngs
RUN make -s
WORKDIR /ngs-tools
RUN ./configure ${BUILD_STYLE} && make -s

FROM alpine:latest
RUN apk add --no-cache libstdc++ libgomp libgcc
COPY --from=build /root/ncbi-outdir/ngs-tools/*/*/*/*/bin /usr/local/ncbi/ngs-tools/bin
ENV PATH=/usr/local/ncbi/ngs-tools/bin:${PATH}
