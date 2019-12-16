#!/usr/bin/env bash

: '
REQUIREMENTS
1. sh interpreter for bash shell. This is already present on GNU/Linux/UNIX machines.
A good alternative on Windows is to use the cmder utility (https://cmder.net/).

2. Python binary on path (named python3, on Windows this can be a link to an executable)
On Windows, create a symbolic link to the Python executable.
Note that python3 can be easily installed on Linux/Unix systems using
an official Python installer or the distro package manager.

3. pipenv installed using pip3.

4. shx installed via the Node Package Manager (ensures cross-platform mkdir and rm)
npm install shx -g

To build:
sh ./run.sh BUILD

To clean:
sh ./run.sh CLEAN
'

if [ "$#" -eq 0 ]; then
    echo "Illegal number of parameters: must be BUILD or CLEAN to run this script"
    echo "USAGE: ./run.sh BUILD to build the figures"
    echo "USAGE: ./run.sh CLEAN to clean the figures"
    exit 1
fi

if [ "$1" = "BUILD" ]
then
    echo "Building all of the figures...this can take time.  Please be patient."
    shx mkdir ../output &&
	  shx mkdir ../output/const-generated/ &&
    shx mkdir ../output/figures-generated &&
	  shx mkdir ../output/tables-generated &&
    pipenv install &&
    pipenv run python3 calibration_figures.py && # DONE
    pipenv run python3 k_q_figure.py &&  # DONE
    pipenv run python3 k_q_example.py &&  # DONE
    pipenv run python3 sp_model_signal_example.py &&  # DONE
    pipenv run python3 dp_model_inv_example.py &&  # DONE
    pipenv run python3 example_signal_processing_figures.py && # DONE
    pipenv run python3 process_all_signal.py &&  # DONE
    pipenv run python3 plots_all_theta_rho.py && # DONE
    pipenv run python3 signal-processing-stats.py &&  # DONE
    pipenv run python3 oat_sensitivity.py
    echo "DONE"
elif [ "$1" = "CLEAN" ]
then
    echo "Cleaning all files created by processing.  Use the BUILD command to re-build the figures."
    pipenv clean &&
	  shx rm -r ../output
    echo "DONE.  Now re-run BUILD to create the figures again."
else
    echo "Option not understood"
fi
