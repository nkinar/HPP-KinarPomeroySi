-----------------------
DATA README
-----------------------
The data in this package is from a paper:

Signal Processing for In-Situ Detection of Effective Heat Pulse Probe Spacing Radius as the Basis of a Self-Calibrating Heat Pulse Probe
by Nicholas J. Kinar, John W. Pomeroy, B.-C. Si

Please see the paper and associated source code for further details.

Data contact:
Dr. Nicholas J. Kinar
Smart Water Systems Lab
University of Saskatchewan
n.kinar@usask.ca

-----------------------

The data in each file is stored in CSV format.  Each file name describes the nominal heat strength used for the heating pulse and the duration of heating.
The total time of experiment is also encoded in the file name.  The associated source code uses the file name to extract the metadata shown below.

YYY-Xwm-Xsec-heating-Xmin-ZZZZ.csv

YYY		is the material type
Xwm 	is the heating strength (W/m)
Xsec	is the time of heating (s)
Xmin	is the total time of experiment (minutes)
ZZZZ	is the repetition identifier

An example of this file format follows:

sand-45wm-8sec-heating-3min-total-rep1.csv

-----------------------

Inside each CSV file, a single time sample is stored per line.  Commas delineate columns.  Column #1
in the table below is the leftmost column in the CSV file.  The time series begins at index 0.
The two sets below show the mapping between column number and variables in each column.

{1, 2, 3, 4, 5, 6, 7, 8, 9}
{Sample number, Vout, dV, I, Rw, T1, T2, DAC tuning code, qprime}

-------------------------------------------------------
Column	|	Identifier
-------------------------------------------------------
1			Time series sample number (starting at 0)
2			Output voltage to the heating element (V)
3			Voltage drop over the current sense resistor connected to the heating element (V)
4			Current through the heating element (A)
5			Resistance of the heating element
6			Temperature inside of the heating needle
7			Temperature of the needle at an offset distance to the heating needle
8			DAC tuning code used for control in the PID loop
9 			Heat strength at each time step

