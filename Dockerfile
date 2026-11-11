# Gen-COMPAS + OpenMM workflow image (GPU at run time; builds without one).
#   docker pull ghcr.io/itaileve-oss/gen-compas:latest
#   docker run --gpus all -it -v $PWD:/workspace ghcr.io/itaileve-oss/gen-compas:latest
#
# Two environments, deliberately separate (torch and OpenMM CUDA pins conflict otherwise):
#   tnfmd       OpenMM 8.4 + PDBFixer + mdtraj   (runs src/gc_md.py, the NAMD stand-in)
#   gen-compas  PyTorch 2.14.1 (cu126) + Gen-COMPAS @ 5004c0d + RiteWeight patch
# Both CUDA user-space stacks come from conda-forge / the torch wheel, so any NVIDIA driver
# with CUDA >= 12.6 works; no nvidia/cuda base image is needed.
FROM mambaorg/micromamba:1.5.10

USER root
RUN apt-get update && apt-get install -y --no-install-recommends git ca-certificates procps \
    && rm -rf /var/lib/apt/lists/*
# workflow configs use /opt/micromamba/envs/<env>/bin/python (the layout on our GPU boxes)
RUN ln -s /opt/conda /opt/micromamba

RUN micromamba create -y -n tnfmd -c conda-forge \
        python=3.12 openmm=8.4.0 cuda-version=12.6 pdbfixer=1.12 mdtraj=1.11.2 \
        numpy scipy pandas matplotlib \
    && micromamba clean -afy

RUN micromamba create -y -n gen-compas -c conda-forge python=3.12 pip \
    && micromamba clean -afy \
    && /opt/conda/envs/gen-compas/bin/pip install --no-cache-dir \
        torch==2.14.1 --index-url https://download.pytorch.org/whl/cu126

ARG GENCOMPAS_COMMIT=5004c0d312c26ee1a8ee0e15ac589228bac5ac2a
COPY patches/riteweight_ergodic_trimming.patch /tmp/riteweight_ergodic_trimming.patch
RUN git clone https://github.com/Tangcyu/Gen-COMPAS.git /opt/Gen-COMPAS \
    && cd /opt/Gen-COMPAS && git checkout ${GENCOMPAS_COMMIT} \
    && git apply /tmp/riteweight_ergodic_trimming.patch \
    && /opt/conda/envs/gen-compas/bin/pip install --no-cache-dir -e . \
    && rm -rf /root/.cache

COPY . /opt/gen-compas-workflow
ENV TORCHDYNAMO_DISABLE=1 \
    PATH=/opt/conda/envs/gen-compas/bin:$PATH

# build-time smoke test (no GPU needed)
RUN /opt/conda/envs/tnfmd/bin/python -c "import openmm, pdbfixer, mdtraj; print('openmm', openmm.__version__)" \
    && /opt/conda/envs/gen-compas/bin/python -c "import torch; print('torch', torch.__version__)" \
    && gen-compas --help > /dev/null && echo "gen-compas CLI OK"

WORKDIR /workspace
CMD ["/bin/bash"]
