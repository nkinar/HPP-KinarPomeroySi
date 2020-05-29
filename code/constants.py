import numpy as np

"""
NOTES:
1. The r0 is stored in ./const-generated as an NPZ file.
The calibrated radius is stored in mm.
"""

# file extensions
PLOT_EXTENSION = '.pdf'
TABLE_EXT = '.txt'
SAVEZ_EXT = '.npz'

# figure size for most square figures
SQUARE_FIGSIZE = (11, 8)
# figure size for most half-figures
HSQUARE_FIGSIZE = (9, 5)
# figure size for long figures
HLONG_FIGSIZE = (14, 5)
HLONG_FIGSIZE_M = (15, 5)
VLONG_FIGSIZE = (17, 12)
LARGE_SQUARE_FIGSIZE = (14, 9)

# font size used for plots
FONT_SIZE = 15

# names of figures
FIG_PATH = '../output/figures-generated/'
TABLE_PATH = '../output/tables-generated/'
CONST_PATH = '../output/const-generated/'
CAL_PLOT_FIGURE_FILENAME = 'cal-plot'
CAL_PLOT_TABLE_FILENAME = 'cal-plot-table'
K_PLOT_TABLE_FILENAME = 'k-plot-table'
K_PLOT_FILENAME = 'k-plot'
CAL_R0_FILENAME = 'cal-r0'
OUTPUT_PROCESSING_FILENAME = 'output-processed'
OUTPUT_PROCESSING_FILENAME_VEC = 'output-processed-vec'
PICKLE_EXT = '.pickle'
K_PLOT_FILENAME_EXAMPLE = 'k-plot-example'
K_PLOT_FILENAME_EXAMPLE_SAND = 'k-plot-example-sand'
SP_SYNTH_FILENAME_EXAMPLE = 'sp-synth-example'
DP_SYNTH_FILENAME_EXAMPLE = 'dp-synth-example'
SAND_EXAMPLE = 'sand-example'
PEAT_EXAMPLE = 'peat-example'
DATA_PICKLE_FILE = CONST_PATH + OUTPUT_PROCESSING_FILENAME + PICKLE_EXT
DATA_PICKLE_FILE_VEC = CONST_PATH + OUTPUT_PROCESSING_FILENAME_VEC + PICKLE_EXT
PLOTS_ALL_THETA_RHO_FILENAME = 'plots-all-theta-rho'
PLOTS_ALL_THERMAL_FILENAME = 'plots-all-k-alpha'
TABLE_NAME_ALL = 'all-comparisons-table'
TABLE_NAME_Q_COMP = 'compare-q-table'

# Constants used to load in the data file
CSV_EXT = '.csv'
CAL_BEGIN = 'cal-'
WM_STR = 'wm-'
SEC_STR = 'sec'
HEATING_STR = '-heating-'
MIN_TOTAL_STR = 'min-total'
DATA_DELIMITER = ','
NO_TEXT = ''
HTML_EXT = '.html'

# Strings related to which curve to downsample or filter
SP_CURVE = 'sp'
DP_CURVE = 'dp'
NO_CURVE = 'none'
BOTH_CURVE = 'both'

# Butterworth filter order
BUTTERWORTH_ORDER = 5       # order of Butterworth filter
BUTTERWORTH_CUTOFF = 10     # [Hz]

# 60 seconds in a minute, used for conversions
SECONDS_IN_MIN = 60.0

# nominal dec places used for output
NOM_DP = 2

# Columns in the data file
NUM_COL = 0
VOUT_COL = 1
dV_COL = 2
I_COL = 3
RW_COL = 4
T1_COL = 5
T2_COL = 6
DAC_CODE_COL = 7
QPRIME_CODE_COL = 8

# step detector window (must be odd)
FWIN_PLAT_WIN = 31
# step detector probability
FWIN_P = 0.01
# decimal places for table formatting
DEC_PLACES_TAB = 2

# Starting value of r0 as set by circuit design
R0_START = 6e-3

# list of q values during calibration
Q_VAL_CAL_LIST = np.asarray([10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75])

# color used for measurements
MEASURED_COLOR = '#f28500'
MODELLED_COLOR = 'limegreen'
PLATEAU_COLOR = 'cornflowerblue'
KNOWN_K_COLOR = 'cornflowerblue'
GREY_COLOR = '#683A5E'
DARK_BLUE = '#0d293e'

# CONSTANTS FOR CALIBRATION
T_COLD_BEGIN = 1.0                          # 1 second at beginning to start temperature
T_HEAT_CAL = 8                              # 8 seconds heating
TIME_TOTAL_CAL = 3 * SECONDS_IN_MIN         # 3 minutes (180 s)
FS_SAMPLE = 120                             # 120 Hz sampling rate
K_WATER = 0.591                             # thermal conductivity of water (W m^-1 K^-1)
RHO_WATER_AGAR = 999                        # density of water with agar gel (kg m^-3)
HEAT_CAPACITY_WATER = 4190                  # heat capacity of water (J kg^-1 K^-1)
DOWNSAMPLE_DP = 1.0                         # downsample dual probe to be 1 s sampling rate
DP_LOWPASS_FILTER_CUT = 10.0                # lowpass filter cut frequency (Hz)
FS_DOWN_DP = 12                             # Downsample for DP

# Volumetric heat capacities of the mineral, organic and water contents (sand)
Cm = 1.9e6          # J m^-3 K^-1
Co = 2.5e6          # J m^-3 K^-1
Cw = 4.186e6        # J m^-3 K^-1

# heat capacities for soil constituents
cm_soil = 0.73e3         # J kg^-3 K^-1
co_soil = 1.9e3          # J kg^-3 K^-1
cw_soil = 4.187e3        # J kg^-3 K^-1

# Volumetric heat capacities for peat
Cm_peat = 1.1e6
Co_peat = 1.0e6

# Densities of the mineral, organic and water contents
rho_m = 2900        # kg m^-3
rho_o = 1300        # kg m^-3
rho_w = 1000        # kg m^-3

# organic content and mineral content for sand
theta_o_sand = 9.2e-3
theta_m_sand = 0.55

# organic content and mineral content for peat
theta_o_peat = 0.30
theta_m_peat = 0.01

QRANGE_IDENT = 'q_range'
R0_IDENT = 'r0'

##################################################

THETA_NOMINAL_SAND = 0.40
DENSITY_NOMINAL_SAND = 1987

THETA_NOMINAL_PEAT = 0.22
DENSITY_NOMINAL_PEAT = 535

# time after the maximum
NOM_TIME_ADD_AFTER_MAX = 2

##################################################
##################################################

SAND_PATH_EXAMPLE = '../hpp-data/hpp-formatted-data/sand/July-21-2017/'
PEAT_PATH_EXAMPLE = '../hpp-data/hpp-formatted-data/peat/July-26-2017/'

##################################################

MAIN_DIR = '../hpp-data/hpp-formatted-data/'
CAL_DIR_NAME = 'calibration'
BIN_DIR_NAME = 'bin'
PEAT_NAME = 'peat'
SAND_NAME = 'sand'
HYPEN_STR = '-'
CSV_STR = '.csv'
SOIL_TYPE_IDENT = 'soil_type'
DATE_IDENT = 'date'

##################################################

# Definitions for vectors
SAND_RADIUS = 'sand-radius'
PEAT_RADIUS = 'peat-radius'
SAND_ORGANIC = 'sand-organic'
PEAT_ORGANIC = 'peat-organic'
SAND_MINERAL = 'sand-mineral'
PEAT_MINERAL = 'peat-mineral'
SAND_TIME = 'sand-time'
PEAT_TIME = 'peat-time'

SAND = 'sand'
PEAT = 'peat'

##################################################

OAT_VEC_FILE = CONST_PATH + 'oat-vec.npz'
SENSITIVITY_ANALYSIS_VEC_FILE = CONST_PATH + 'sensitivity-analysis.pickle'

OAT_SENSITIVITY_FIG_NAME_THETA = FIG_PATH + 'oat-theta' + PLOT_EXTENSION
OAT_SENSITIVITY_FIG_NAME_RHO = FIG_PATH + 'oat-rho' + PLOT_EXTENSION

##################################################




