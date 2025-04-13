

Basic Install Instructions

1) Download Prob3plusplus:
https://github.com/rogerwendell/Prob3plusplus

2) Install Prob3plusplus (simply go inside and type 'make')

3) Install pyROOT.
	- $ conda config --set channel_priority strict
	- $ conda create -c conda-forge --name <my-environment> root
	- $ conda activate <my-environment>

	Also install in this environment: fire, matplotlib, cython
	-use: conda install packagename -c conda-forge

4) Build cython wrapper (this assumes you are working on a macOS machine with Apple Silicon arch).
	- Edit hardcoded paths within wrapper_Prob3 "\*.pyx" files.
	- Edit hardcoded paths within wrapper_Prob3 setup.py file.
	- Run the command: python setup.py build    (build is the name of the folder where the libraries will be stored)
    - If you encounter miss-match versions during the linking process, you might need to use: export MACOSX_DEPLOYMENT_TARGET=11.0 (change version as needed).
	- Copy .so files from inside wrapper_Prob3/build/lib into /build

5) Change hardcoded paths within configuration files

6) Create some setup script that activates your conda environment and adds the project to your python path: export PYTHONPATH=$PYTHONPATH:/path/to/project/folder/
