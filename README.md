Step 1:

```sudo docker build -t quick-surge .```

Step 2: 

```sudo docker run -v   /mnt/DATA/QuickSurge-docker/QuickSurge:/QuickSurge -ti --rm quick-surge /bin/bash -c ". /opt/miniconda/bin/activate; cd /QuickSurge/Python_Codes; python main_IDA_PDNA_StormSurge.py"```
