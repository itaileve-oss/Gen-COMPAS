"""Diffusion diagnostics (control folder only). Train a fresh TNF diffusion model with
the iteration-0 effective config, changing ONLY the listed settings, then sample.
Usage: exp_train.py OUT_DIR EPOCHS DCD [noise1,noise2,...]"""
import sys, copy, os, yaml
sys.path.insert(0, '/workspace/gc/external/Gen-COMPAS')
from common.diffusion_train import train_diffusion_model
from common.diffusion_sample import run_diffusion_inference
out, epochs, dcd = sys.argv[1], int(sys.argv[2]), sys.argv[3]
noises = [float(x) for x in (sys.argv[4] if len(sys.argv) > 4 else "10,2").split(",")]
os.makedirs(out, exist_ok=True)
cfg = yaml.safe_load(open('/workspace/gc/Iterations/0th.Iteration/effective_config.yaml'))
g = copy.deepcopy(cfg['Generative'])
g['save_dir'] = os.path.join(out, 'models'); g['init_checkpoint_path'] = None
g['training']['epochs'] = epochs; g['training']['save_interval'] = max(epochs, 1)
g['data']['dcd_path'] = dcd
print('TRAIN', epochs, 'epochs on', dcd, '->', g['save_dir'], flush=True)
train_diffusion_model(g)
for nz in noises:
    s = copy.deepcopy(g)
    s['inference']['checkpoint'] = os.path.join(g['save_dir'], 'best_model.pt')
    s['inference']['noise_scale'] = nz; s['inference']['num_samples'] = 1000
    s['inference']['output'] = os.path.join(out, f'candidates_noise{nz:g}.dcd')
    print('SAMPLE noise', nz, flush=True)
    run_diffusion_inference(s)
