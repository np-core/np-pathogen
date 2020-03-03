Bootstrap: docker
From: continuumio/miniconda3

%labels
    name np-sepsis
    version v0.1.0
    author esteinig

%post
    export PATH=/opt/conda/bin:$PATH

    apt-get update && apt install libgl1-mesa-glx qt5-default build-essential -y

    conda install -c conda-forge -c bioconda --yes \
        trimmomatic \
        pysam \
        pilon=1.23 \
        racon \
        minimap2 \
        shovill \
        flye \
        nanofilt \
        nanostat \
        bandage \
        kraken2 \
        pip \
        pyqt \
        && pip install /nanopath \
        && conda clean -a \
        && find /opt/conda/ -follow -type f -name '*.a' -delete \
        && find /opt/conda/ -follow -type f -name '*.pyc' -delete

    sed -i '16s/512m/8g/' /opt/conda/share/pilon-1.23-2/pilon
    sed -i '16s/1g/16g/' /opt/conda/share/pilon-1.23-2/pilon
%files
    ../nanopath /nanopath

%environment
    export PATH=:/opt/conda/bin:$PATH