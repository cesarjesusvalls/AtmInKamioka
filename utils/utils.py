import os, glob, yaml
def list_available_devices():
    devs = dict(cpu=torch.device('cpu'))

    if torch.cuda.is_available():
        devs['cuda'] = torch.device('cuda:0')

    if torch.backends.mps.is_available():
        devs['mps'] = torch.device('mps')

    return devs
def get_device(request):

    devs = list_available_devices()

    if not request in devs:
        print(request,'not supported')
        return None
    else:
        return devs[request]

def get_config_dir():
    return os.path.join(os.path.dirname(__file__),'../config')

def list_config(full_path=False):

    fs = glob.glob(os.path.join(get_config_dir(), '*.yaml'))

    if full_path:
        return fs

    return [os.path.basename(f)[:-5] for f in fs]

def get_config(name):

    options = list_config()
    results = list_config(True)

    if name in options:
        return results[options.index(name)]

    alt_name = name + '.yaml'
    if alt_name in options:
        return results[options.index(alt_name)]

    print('No data found for config name:',name)
    raise NotImplementedError
def load_config(name:str):

    return yaml.safe_load(open(get_config(name),'r'))