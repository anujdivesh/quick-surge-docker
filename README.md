sudo docker build -it quick-surge .
sudo docker run -v /mnt/DATA/QuickSurge-docker/QuickSurge:/QuickSurge -ti --rm quick-surge /bin/bash -c ". /opt/miniconda/bin/activate; cd /QuickSurge/Python_Codes; python main_IDA_PDNA_StormSurge.py"
