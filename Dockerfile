FROM ubuntu:22.04

RUN DEBIAN_FRONTEND=noninteractive apt-get -y update 
ENV TZ=Etc/UTC
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get -y install build-essential g++ gfortran mpich gcc python3-dev python3-pip virtualenv wget zlib1g-dev vim htop libnetcdf-dev libnetcdff-dev bzip2 cmake cpio curl git gosu libblas-dev liblapack-dev libmpich-dev 
#RUN apt-get -y install apache2

RUN mkdir -p /opt/conda 

RUN mkdir -p /QuickSurge
COPY ./ /QuickSurge/
WORKDIR /QuickSurge 
#VOLUME /model

RUN bash /QuickSurge/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh -b -p /opt/miniconda 

RUN . /opt/miniconda/bin/activate && conda env update --name base --file /QuickSurge/environment.yml 

#listen to port 80
#EXPOSE 80

#start Apache
#CMD service apache2 start && \
#    tail -F -n0 /etc/hosts
#EXPOSE 5003
#CMD ["/opt/miniconda/bin/conda", "run", "--no-capture-output", "-n", "base", "python" "/model/Python_Codes/main_IDA_PDNA_StormSurge.py"]
#ENTRYPOINT  [". /opt/miniconda/bin/activate", "cd /model/Python_Codes", "python main_IDA_PDNA_StormSurge.py"]
CMD [ "/bin/bash" ]