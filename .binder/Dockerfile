FROM rocker/binder:4.3.0
LABEL maintainer='viv3kanand'
COPY --chown=${NB_USER} . ${HOME}
USER ${NB_USER}



RUN wget https://github.com/viv3kanand/SCA12-RNA-Seq-Analysis/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2023-06-09/')); devtools::install_deps()"

RUN rm DESCRIPTION.1; exit 0
