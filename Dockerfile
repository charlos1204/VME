FROM rocker/shiny:4.0.3

LABEL maintainer="Carlos Armando Garcia Perez <carlos.garcia@helmholtz-muenchen.de>"
LABEL version="1.0"
LABEL description="Shiny app of the Community Explorer"

# Install git, wget, python-dev, pip, BLAS + LAPACK and other dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
  nano \
  g++ \
  wget \
  libxml2-dev \
  software-properties-common

RUN apt-get update && apt-get install -y \
  libcurl4-openssl-dev libglpk-dev && \
  mkdir -p /var/lib/shiny-server/bookmarks/shiny

COPY installapps.R .

# install libraries
RUN Rscript installapps.R

ADD datasets /root/datasets
COPY app /root/app
COPY Rprofile.site /usr/local/lib/R/etc/Rprofile.site

RUN chmod -R 755 /root/app
RUN chmod -R 755 /usr/local/lib/R/etc

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/app')"]

#sudo docker image ls
#sudo docker system prune -a
#sudo docker build -t charlos1204/vme .
#sudo docker push charlos1204/vme
#sudo docker run -p 3838:3838 charlos1204/vme

