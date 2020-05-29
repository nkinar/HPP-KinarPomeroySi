# Signal Processing for In-Situ Detection of Effective Heat Pulse Probe Spacing Radius as the Basis of a Self-Calibrating Heat Pulse Probe
## by Nicholas J. Kinar, John W. Pomeroy, B.-C. Si

*Journal of Geoscientific Instrumentation, Methods and Data Systems (GI)*

## What is it?

This is the complete Python source code required to reproduce the paper.  The source code also produces data figures and comparisons.

## Paper Abstract
A sensor comprised of an electronic circuit and a hybrid single and dual heat pulse probe was constructed and tested along with a novel signal processing procedure to determine changes in the effective dual-probe spacing radius over the time of measurement.  The circuit utilized a proportional–integral–derivative (PID) controller to control heat inputs into the soil medium in lieu of a variable resistor.  The system was designed for on-board signal processing and implemented USB, RS-232 and SDI-12 interfaces for Machine-to-Machine (M2M) exchange of data, thereby enabling heat inputs to be adjusted to soil conditions and data availability shortly after the time of experiment.  Signal processing was introduced to provide a simplified single-probe model to determine thermal conductivity instead of reliance on late-time logarithmic curve-fitting.  Homomorphic and derivative filters were used with a dual-probe model to detect changes in the effective probe spacing radius over the time of experiment to compensate for physical changes in radius as well as model and experimental error.  Theoretical constraints were developed for an efficient inverse of the exponential integral on an embedded system.  Application of the signal processing to experiments on sand and peat improved the estimates of soil water content and bulk density compared to methods of curve-fitting nominally used for heat pulse probe experiments.  Applications of the technology may be especially useful for soil and environmental conditions where effective changes in probe spacing radius need to be detected and compensated over the time of experiment.

## Code Dependencies

* Python 3
* numpy
* scipy
* matplotlib
* prettytable
* natsort
* beautifulsoup4
* pandas
* lxml

For the `run.sh` script:
* pipenv
* `shx` installed via the Node Package Manager as `npm install shx -g`
* ImageMagick for figure conversions

All of the Python dependencies can be installed using `pip3` if the `run.sh` script is not to be used.
When using `pipenv`, the dependencies will be automatically installed by the `run.sh` script.

An `sh` interpreter is also required for the shell script `run.sh`. This is already present on `GNU/Linux/UNIX` machines
including `Mac OS X`.  A good alternative for Windows is to use the `cmder` utility
downloaded from <https://cmder.net/>.

It is nominally assumed that the interpreter is `bash`.

## Use

The shell script `run.sh` will automatically create a virtual environment and
install all required packages using the excellent `pipenv`.

To build all files:

```
sh ./run.sh BUILD
```
To clean all processing and to delete generated files:

```
sh ./run.sh CLEAN
```

## Data DOI

The data and code is available in this repo.  See `readme.txt` in the `hpp-data` folder.

The data can also be downloaded from Figshare as <y>.

## Built With

* Python 3
* pipenv

## Built @

* Smart Water Systems Lab (University of Saskatchewan)

## Author

* **Nicholas J. Kinar** - code and mathematics

## License

This project is licensed under the GPLv3 License - see the [LICENSE.md](LICENSE.md) file for details.
The GPLv3 License ensures that the code can be freely shared.
