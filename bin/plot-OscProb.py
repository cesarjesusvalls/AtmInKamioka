import yaml
import os
import fire
from app import oscProb
from utils import list_config, load_config

def main(config_file):
    cfg = dict()
    if not os.path.isfile(config_file) and config_file in list_config():
        cfg=load_config(config_file)
    else:
        with open(config_file,'r') as f:
            cfg=yaml.safe_load(f)

    list_config('config')
    pltOscProb = oscProb(cfg)
    pltOscProb._load_barger_prop()

    pltOscProb._plot_osc_prob()

if __name__ == '__main__':
    fire.Fire(main)