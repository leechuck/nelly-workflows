# Self-contained Exomiser image: JRE + the exomiser-cli 15.1.0 distribution + a
# real wrapper script. The official exomiser/exomiser-cli image is distroless
# (no shell, fixed jib entrypoint), which can't run our setup script -- so we
# build our own from a JRE base that has bash.
#
# Since v4 the distribution zip is downloaded from the GitHub release at build
# time (asset verified: exomiser-cli-15.1.0-distribution.zip, ~170 MB, top-level
# dir exomiser-cli-15.1.0/), instead of COPYing a pre-unpacked copy from the
# build context.
FROM eclipse-temurin:21-jre

RUN apt-get update && apt-get install -y --no-install-recommends \
        curl unzip gzip grep coreutils && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /exomiser && \
    curl -fsSL -o /tmp/exomiser-cli-15.1.0-distribution.zip \
      https://github.com/exomiser/Exomiser/releases/download/15.1.0/exomiser-cli-15.1.0-distribution.zip && \
    unzip -q /tmp/exomiser-cli-15.1.0-distribution.zip -d /tmp/ex && \
    mv /tmp/ex/exomiser-cli-15.1.0/* /exomiser/ && \
    rm -rf /tmp/ex /tmp/exomiser-cli-15.1.0-distribution.zip

COPY run-exomiser.sh /run-exomiser.sh
RUN chmod +x /run-exomiser.sh

WORKDIR /work
# No ENTRYPOINT: CWL supplies the full command (baseCommand: [bash,
# /run-exomiser.sh]). An ENTRYPOINT would be PREPENDED to the CWL command by the
# Arvados runtime, doubling "bash /run-exomiser.sh" and shifting the positional
# args (the data path landed in $3/assembly).
CMD ["bash"]
