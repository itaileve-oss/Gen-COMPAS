"""Resample the trained TNF v2 iteration-0 diffusion model with another noise scale.
Reads the iteration's effective config + model READ-ONLY; writes only to OUT_DIR."""
import sys, copy, os, yaml
sys.path.insert(0, '/workspace/gc/external/Gen-COMPAS')
from common.diffusion_sample import run_diffusion_inference
noise, out_dir = float(sys.argv[1]), sys.argv[2]
os.makedirs(out_dir, exist_ok=True)
cfg = yaml.safe_load(open('/workspace/gc/Iterations/0th.Iteration/effective_config.yaml'))
g = copy.deepcopy(cfg['Generative'])
g['inference']['noise_scale'] = noise
g['inference']['num_samples'] = 1000
g['inference']['output'] = os.path.join(out_dir, f'candidates_noise{noise:g}.dcd')
print('checkpoint', g['inference'].get('checkpoint'), '| save_dir', g.get('save_dir'), '| out', g['inference']['output'])
run_diffusion_inference(g)
