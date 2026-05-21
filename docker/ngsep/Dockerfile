FROM eclipse-temurin:21

RUN wget https://github.com/NGSEP/NGSEPcore/releases/download/v5.1.0/NGSEPcore_5.1.0.jar -O /usr/local/bin/NGSEPcore.jar

ENTRYPOINT [ "java", "-jar", "/usr/local/bin/NGSEPcore.jar" ]