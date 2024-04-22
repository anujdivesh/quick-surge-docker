#!/bin/bash

docker run -v /mnt/DATA/QuickSurge-docker/QuickSurge:/QuickSurge --rm quick-surge /bin/bash -c ". /opt/miniconda/bin/activate; cd /QuickSurge/Python_Codes; python main_IDA_PDNA_StormSurge.py"
echo "Model run Successful!"
