#!/usr/bin/env bash
# Provision a bare vast.ai GPU box for Gen-COMPAS on TNF. Run ON the instance.
# Two envs, deliberately separate (torch and openmm CUDA pins fight otherwise):
#   tnfmd       OpenMM + pdbfixer + mdtraj  (MD; gc_md.py runs here)
#   gen-compas  torch(CUDA) + Gen-COMPAS    (diffusion, VCN, RiteWeight)
set -euo pipefail
GC=/workspace/gc
MAMBA_ROOT=/opt/micromamba
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
DRIVER_CUDA=$(nvidia-smi | sed -n 's/.*CUDA Version: \([0-9.]*\).*/\1/p' | head -1)
echo "driver CUDA ceiling ${DRIVER_CUDA}"
if [ ! -x "${MAMBA_ROOT}/bin/micromamba" ]; then
  mkdir -p "${MAMBA_ROOT}/bin"
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "${MAMBA_ROOT}" bin/micromamba
fi
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT}"
eval "$("${MAMBA_ROOT}/bin/micromamba" shell hook -s bash)"
# ALWAYS pin cuda-version: an unpinned openmm-family install once pulled 13.5 above a
# 13.2 driver and PME silently fell back to CPU while a toy smoke test still passed.
micromamba create -y -n tnfmd -c conda-forge python=3.12 openmm "cuda-version=12.6" pdbfixer mdtraj numpy
micromamba create -y -n gen-compas -c conda-forge python=3.12 libstdcxx-ng libgcc-ng pip git
micromamba run -n gen-compas pip install -q torch --index-url https://download.pytorch.org/whl/cu126
cd "${GC}/external" && [ -d Gen-COMPAS ] || git clone -q https://github.com/Tangcyu/Gen-COMPAS
cd "${GC}/external/Gen-COMPAS" && git checkout -q 5004c0d && micromamba run -n gen-compas pip install -q -e .
micromamba run -n gen-compas python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"
micromamba run -n gen-compas gen-compas --help | head -3
echo "=== ready: next run gencompas/build_endpoints.sh (real PME benchmark = prep ns/day; also watch nvidia-smi) ==="
