import numpy as np
import pickle
from gen_mid_vec import gen_mid_vec
from process_all_signal import process_all, create_theta_rho_data_vec
from float_round import float_round_str
from constants import *
from comparisons import compute_rmse, compute_mb
from AutoDict import AutoDict
from constant_labels import create_label

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('ggplot')
mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams.update({'font.size': 12})


def get_oat_filename_processing(val, add):
    """
    Get the filename required for initial processing
    :return:
    """
    fn = CONST_PATH + OUTPUT_PROCESSING_FILENAME + float_round_str(val, 4) + add + PICKLE_EXT
    return fn


def get_oat_filename_processing_vec(val, add):
    """
    Get the filename required for vector processing
    :param val:
    :return:
    """
    fn = CONST_PATH + OUTPUT_PROCESSING_FILENAME_VEC + float_round_str(val, 4) + add + PICKLE_EXT
    return fn


def run_all_process(check_if_file_exists=False):
    """
    Run the processing for the OAT sensitivity
    :return:
    """
    # number of points over which to conduct the sensitivity analysis
    n = 10

    # change in delta r as the starting value of r
    rmid = R0_START
    rmin = 5e-3
    rmax = 14e-3
    delta_r = gen_mid_vec(rmin, rmid, rmax, n)

    print('Processing for sand radius...')
    for sand_radius_test in delta_r:
        print('r = ', sand_radius_test)
        df0 = get_oat_filename_processing(sand_radius_test, SAND_RADIUS)
        process_all(df0,
                    check_if_file_exists=check_if_file_exists,
                    r_nominal=sand_radius_test,
                    theta_o=None,
                    theta_m=None,
                    add_time_after_max=None, sand=True, peat=False)
        df1 = get_oat_filename_processing_vec(sand_radius_test, SAND_RADIUS)
        create_theta_rho_data_vec(df0, df1, check_if_file_exists)

    print('Processing for the peat radius...')
    for peat_radius_test in delta_r:
        print('r = ', peat_radius_test)
        df0 = get_oat_filename_processing(peat_radius_test, PEAT_RADIUS)
        process_all(df0,
                    check_if_file_exists=check_if_file_exists,
                    r_nominal=peat_radius_test,
                    theta_o=None,
                    theta_m=None,
                    add_time_after_max=None, sand=False, peat=True)
        df1 = get_oat_filename_processing_vec(peat_radius_test, PEAT_RADIUS)
        create_theta_rho_data_vec(df0, df1, check_if_file_exists)

    ########################################################
    # ORGANIC CONTENT
    ########################################################

    # [1]
    # theta_o for sand organic content
    print('Processing for sand organic content...')
    theta_o_min_s = 1.0e-3
    theta_o_mid_s = theta_o_sand    # 9.2e-3
    theta_o_high_s = 0.2
    delta_theta_o_sand = gen_mid_vec(theta_o_min_s, theta_o_mid_s, theta_o_high_s, n)
    for theta_o_sand_test in delta_theta_o_sand:
        print('theta_o = ', theta_o_sand_test)
        df0 = get_oat_filename_processing(theta_o_sand_test, SAND_ORGANIC)
        process_all(df0,
                    check_if_file_exists=check_if_file_exists,
                    r_nominal=None,
                    theta_o=theta_o_sand_test,
                    theta_m=None,
                    add_time_after_max=None,
                    sand=True, peat=False)
        df1 = get_oat_filename_processing_vec(theta_o_sand_test, SAND_ORGANIC)
        create_theta_rho_data_vec(df0, df1, check_if_file_exists)

    # [2]
    # theta_o for the organic content for peat
    print('Processing for peat organic content...')
    theta_o_min_p = 0.1
    theta_o_mid_p = theta_o_peat    # 0.49
    theta_o_high_p = 0.80
    delta_theta_o_peat = gen_mid_vec(theta_o_min_p, theta_o_mid_p, theta_o_high_p, n)
    for theta_o_peat_test in delta_theta_o_peat:
        print('theta_o = ', theta_o_peat_test)
        df0 = get_oat_filename_processing(theta_o_peat_test, PEAT_ORGANIC)
        process_all(df0,
                    check_if_file_exists=check_if_file_exists,
                    r_nominal=None,
                    theta_o=theta_o_peat_test,
                    theta_m=None,
                    add_time_after_max=None,
                    sand=False, peat=True)
        df1 = get_oat_filename_processing_vec(theta_o_peat_test, PEAT_ORGANIC)
        create_theta_rho_data_vec(df0, df1, check_if_file_exists)

    ########################################################
    # MINERAL CONTENT
    ########################################################

    # [3]
    # theta_m for sand
    print('Processing for sand mineral content...')
    theta_m_min_s = 0.3
    theta_m_mid_s = theta_m_sand    # 0.59
    theta_m_high_s = 0.80
    delta_theta_m_sand = gen_mid_vec(theta_m_min_s, theta_m_mid_s, theta_m_high_s, n)
    for theta_m_sand_test in delta_theta_m_sand:
        print('theta_m = ', theta_m_sand_test)
        df0 = get_oat_filename_processing(theta_m_sand_test, SAND_MINERAL)
        process_all(df0,
                    check_if_file_exists=check_if_file_exists,
                    r_nominal=None,
                    theta_o=None,
                    theta_m=theta_m_sand_test,
                    add_time_after_max=None,
                    sand=True, peat=False)
        df1 = get_oat_filename_processing_vec(theta_m_sand_test, SAND_MINERAL)
        create_theta_rho_data_vec(df0, df1, check_if_file_exists)

    # [4]
    # theta_m for peat
    print('Processing for peat mineral content...')
    theta_m_min_p = 0.001
    theta_m_mid_p = theta_m_peat    # 0.01
    theta_m_high_p = 0.2
    delta_theta_m_peat = gen_mid_vec(theta_m_min_p, theta_m_mid_p, theta_m_high_p, n)
    for theta_m_peat_test in delta_theta_m_peat:
        print('theta_m = ', theta_m_peat_test)
        df0 = get_oat_filename_processing(theta_m_peat_test, PEAT_MINERAL)
        process_all(df0,
                    check_if_file_exists=check_if_file_exists,
                    r_nominal=None,
                    theta_o=None,
                    theta_m=theta_m_peat_test,
                    add_time_after_max=None,
                    sand=False, peat=True)
        df1 = get_oat_filename_processing_vec(theta_m_peat_test, PEAT_MINERAL)
        create_theta_rho_data_vec(df0, df1, check_if_file_exists)

    #####################################################################################
    # TIME AFTER / BEFORE MAX (can be negative since this is relative to the peak time)
    #####################################################################################
    # [5]
    low_time = -6.0
    mid_time = NOM_TIME_ADD_AFTER_MAX
    high_time = 6.0
    delta_time = gen_mid_vec(low_time, mid_time, high_time, n)

    print('Processing time delay for sand...')
    for time_sand_test in delta_time:
        print('time delay = ', time_sand_test)
        df0 = get_oat_filename_processing(time_sand_test, SAND_TIME)
        process_all(df0,
                    check_if_file_exists=check_if_file_exists,
                    r_nominal=None,
                    theta_o=None,
                    theta_m=None,
                    add_time_after_max=time_sand_test,
                    sand=True, peat=False)
        df1 = get_oat_filename_processing_vec(time_sand_test, SAND_TIME)
        create_theta_rho_data_vec(df0, df1, check_if_file_exists)

    print('Processing time for peat...')
    for time_sand_test in delta_time:
        print('time delay = ', time_sand_test)
        df0 = get_oat_filename_processing(time_sand_test, PEAT_TIME)
        process_all(df0,
                    check_if_file_exists=check_if_file_exists,
                    r_nominal=None,
                    theta_o=None,
                    theta_m=None,
                    add_time_after_max=time_sand_test,
                    sand=False, peat=True)
        df1 = get_oat_filename_processing_vec(time_sand_test, PEAT_TIME)
        create_theta_rho_data_vec(df0, df1, check_if_file_exists)

        print('Saving out OAT vector file...')
        np.savez(OAT_VEC_FILE, delta_r=delta_r,
                 delta_theta_o_sand=delta_theta_o_sand,
                 delta_theta_o_peat=delta_theta_o_peat,
                 delta_theta_m_sand=delta_theta_m_sand,
                 delta_theta_m_peat=delta_theta_m_peat,
                 delta_time=delta_time)


def load_file_compute_rmsd_mb(fn):
    """
    Load the file and compute the required outputs
    :param fn:
    :param theta_comp_sand:
    :param rho_comp_sand:
    :param theta_comp_peat:
    :param rho_comp_peat:
    :return:
    """
    experiments = pickle.load(open(fn, "rb"))
    name = SAND_NAME
    if SAND_NAME in experiments:
        name = SAND_NAME
    elif PEAT_NAME in experiments:
        name = PEAT_NAME
    print('len = ', len(experiments[name]))
    theta_nom, rho_nom, theta_sig, rho_sig, \
    alpha_det, alpha_sig, k_det, kdet_sig, kdet_linear, _qav, _qknown = experiments[name]
    return theta_nom, rho_nom, theta_sig, rho_sig, \
    alpha_det, alpha_sig, k_det, kdet_sig, kdet_linear  # note that the qav and qknown are not returned here


def load_each_assemble(val, ident_str):
    """
    Load each assembly
    :param val:
    :param ident_str:
    :return:
    """
    fn = get_oat_filename_processing_vec(val, ident_str)
    print('fn = ', fn)
    theta_nom, rho_nom, theta_sig, rho_sig, \
    alpha_det, alpha_sig, k_det, kdet_sig, kdet_linear = load_file_compute_rmsd_mb(fn)
    return theta_nom, rho_nom, theta_sig, rho_sig, \
    alpha_det, alpha_sig, k_det, kdet_sig, kdet_linear


def compute_over_vec(delta_vec, ident_vec):

    dictionary = AutoDict()

    for idx, ident_str in enumerate(ident_vec):                     # iterate over identify vector

        # identification string and the index
        print('ident = ', ident_str, 'idx = ', idx)

        if SAND in ident_str:
            theta_comp = THETA_NOMINAL_SAND
            rho_comp = DENSITY_NOMINAL_SAND
        elif PEAT in ident_str:
            theta_comp = THETA_NOMINAL_PEAT
            rho_comp = DENSITY_NOMINAL_PEAT
        else:
            raise ValueError('Wrong identification string')

        # vector with delta change
        dvec = delta_vec[idx]

        # create vectors
        theta_nom_rmsd_vec = []
        theta_nom_mb_vec = []
        rho_nom_rmse_vec = []
        rho_nom_mb_vec = []

        theta_sig_rmsd_vec = []
        theta_sig_mb_vec = []
        rho_sig_rmse_vec = []
        rho_sig_mb_vec = []

        for d in dvec:                                              # iterate over delta vector value

            # print the value over which the vector has occurred
            print('d = ', d)

            theta_nom, rho_nom, theta_sig, rho_sig, \
            alpha_det, alpha_sig, k_det, kdet_sig, kdet_linear = load_each_assemble(d, ident_str)

            # compute the RMSD and MB
            theta_nom_rmsd = compute_rmse(theta_comp, theta_nom)
            theta_nom_mb = compute_mb(theta_comp, theta_nom)
            rho_nom_rmse = compute_rmse(rho_comp, rho_nom)
            rho_nom_mb = compute_mb(rho_comp, rho_nom)

            theta_sig_rmsd = compute_rmse(theta_comp, theta_sig)
            theta_sig_mb = compute_mb(theta_comp, theta_sig)
            rho_sig_rmse = compute_rmse(rho_comp, rho_sig)
            rho_sig_mb = compute_mb(rho_comp, rho_sig)

            theta_nom_rmsd_vec.append(theta_nom_rmsd)
            theta_nom_mb_vec.append(theta_nom_mb)
            rho_nom_rmse_vec.append(rho_nom_rmse)
            rho_nom_mb_vec.append(rho_nom_mb)

            theta_sig_rmsd_vec.append(theta_sig_rmsd)
            theta_sig_mb_vec.append(theta_sig_mb)
            rho_sig_rmse_vec.append(rho_sig_rmse)
            rho_sig_mb_vec.append(rho_sig_mb)

        # store in a dictionary
        # dvec = change in quantity (such as r)
        dictionary[ident_str] = (dvec, theta_nom_rmsd_vec,  theta_nom_mb_vec, rho_nom_rmse_vec, rho_nom_mb_vec,
                                 theta_sig_rmsd_vec, theta_sig_mb_vec, rho_sig_rmse_vec, rho_sig_mb_vec)
    # DONE loop
    return dictionary
# DONE


def compute_rmsd_mb_oat():
    """
    Compute the RMSD and MB of all the data (run the functions above)
    :return:
    """
    print('Running analysis to compute the RMSD and the MB...')
    vec = np.load(OAT_VEC_FILE)
    delta_r = vec['delta_r']
    delta_theta_o_sand = vec['delta_theta_o_sand']
    delta_theta_o_peat = vec['delta_theta_o_peat']
    delta_theta_m_sand = vec['delta_theta_m_sand']
    delta_theta_m_peat = vec['delta_theta_m_peat']
    delta_time = vec['delta_time']

    # NOTE that both the IDENT_VEC and the NUM_VEC must coincide
    # NOTE that the delta_r at the beginning and the delta_time at the end must be repeated twice since these
    # vectors are used twice as well
    NUM_VEC = [delta_r, delta_r, delta_theta_o_sand, delta_theta_o_peat, delta_theta_m_sand, delta_theta_m_peat,
               delta_time, delta_time]
    IDENT_VEC = [SAND_RADIUS, PEAT_RADIUS, SAND_ORGANIC, PEAT_ORGANIC, SAND_MINERAL, PEAT_MINERAL, SAND_TIME, PEAT_TIME]

    dictionary = compute_over_vec(NUM_VEC, IDENT_VEC)
    print('Saving dictionary out to pickle file...')
    pickle.dump(dictionary, open(SENSITIVITY_ANALYSIS_VEC_FILE, 'wb'))
    print('DONE saving to pickle file')


def plot_oat(show=False):
    """
    Plot for sensitivity analysis
    :return:
    """
    # load in the data from the dictionary file
    print('Loading in the data...')
    dictionary = pickle.load(open(SENSITIVITY_ANALYSIS_VEC_FILE, "rb"))
    print('DONE loading in the data.')

    #######################################################################################

    # RADIUS
    rvec_sand, theta_nom_rmsd_vec_sand_radius, theta_nom_mb_vec_sand_radius, rho_nom_rmse_vec_sand_radius, \
    rho_nom_mb_vec_sand_radius, theta_sig_rmsd_vec_sand_radius, theta_sig_mb_vec_sand_radius, \
    rho_sig_rmse_vec_sand_radius, rho_sig_mb_vec_sand_radius = dictionary[SAND_RADIUS]

    rvec_peat, theta_nom_rmsd_vec_peat_radius, theta_nom_mb_vec_peat_radius, rho_nom_rmse_vec_peat_radius, \
    rho_nom_mb_vec_peat_radius, theta_sig_rmsd_vec_peat_radius, theta_sig_mb_vec_peat_radius, \
    rho_sig_rmse_vec_peat_radius, rho_sig_mb_vec_peat_radius = dictionary[PEAT_RADIUS]

    # ORGANIC CONTENT
    dvec_sand_theta_o, theta_nom_rmsd_vec_sand_theta_o, theta_nom_mb_vec_sand_theta_o, rho_nom_rmse_vec_sand_theta_o, \
    rho_nom_mb_vec_sand_theta_o, theta_sig_rmsd_vec_sand_theta_o, theta_sig_mb_vec_sand_theta_o, \
    rho_sig_rmse_vec_sand_theta_o, rho_sig_mb_vec_sand_theta_o = dictionary[SAND_ORGANIC]

    dvec_peat_theta_o, theta_nom_rmsd_vec_peat_theta_o, theta_nom_mb_vec_peat_theta_o, rho_nom_rmse_vec_peat_theta_o, \
    rho_nom_mb_vec_peat_theta_o, theta_sig_rmsd_vec_peat_theta_o, theta_sig_mb_vec_peat_theta_o, \
    rho_sig_rmse_vec_peat_theta_o, rho_sig_mb_vec_peat_theta_o = dictionary[PEAT_ORGANIC]

    # MINERAL CONTENT
    (dvec_sand_theta_m, theta_nom_rmsd_vec_sand_theta_m, theta_nom_mb_vec_sand_theta_m, rho_nom_rmse_vec_sand_theta_m, rho_nom_mb_vec_sand_theta_m,
     theta_sig_rmsd_vec_sand_theta_m, theta_sig_mb_vec_sand_theta_m, rho_sig_rmse_vec_sand_theta_m,
     rho_sig_mb_vec_sand_theta_m) = dictionary[SAND_MINERAL]

    (dvec_peat_theta_m, theta_nom_rmsd_vec_peat_theta_m, theta_nom_mb_vec_peat_theta_m, rho_nom_rmse_vec_peat_theta_m, rho_nom_mb_vec_peat_theta_m,
     theta_sig_rmsd_vec_peat_theta_m, theta_sig_mb_vec_peat_theta_m, rho_sig_rmse_vec_peat_theta_m,
     rho_sig_mb_vec_peat_theta_m) = dictionary[PEAT_MINERAL]

    # TIME SHIFT
    (dvec_sand_time, theta_nom_rmsd_vec_sand_time, theta_nom_mb_vec_sand_time,
     rho_nom_rmse_vec_sand_time, rho_nom_mb_vec_sand_time,
     theta_sig_rmsd_vec_sand_time, theta_sig_mb_vec_sand_time,
     rho_sig_rmse_vec_sand_time, rho_sig_mb_vec_sand_time) = dictionary[SAND_TIME]

    (dvec_peat_time, theta_nom_rmsd_vec_peat_time, theta_nom_mb_vec_peat_time,
     rho_nom_rmse_vec_peat_time, rho_nom_mb_vec_peat_time,
     theta_sig_rmsd_vec_peat_time, theta_sig_mb_vec_peat_time,
     rho_sig_rmse_vec_peat_time, rho_sig_mb_vec_peat_time) = dictionary[PEAT_TIME]

    #######################################################################################

    # There will be two plots:
    # plot theta (3 x 2 = 6)
    # plot rho (3 x 2 = 6)

    #########################
    # THETA ONLY
    #########################

    fig = plt.figure(num=None, figsize=VLONG_FIGSIZE)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # (a)
    ax = fig.add_subplot(4, 2, 1)
    ax.plot(rvec_sand, theta_nom_rmsd_vec_sand_radius, label='Heating and Cooling DP Sand')
    ax.plot(rvec_sand, theta_sig_rmsd_vec_sand_radius, label='Signal Processing SP and DP Sand')
    ax.plot(rvec_peat, theta_nom_rmsd_vec_peat_radius, label='Heating and Cooling DP Peat')
    ax.plot(rvec_peat, theta_sig_rmsd_vec_peat_radius, label='Signal Processing SP and DP Peat')
    ax.set_xlabel(create_label('$r$', 'mm'))
    ax.set_ylabel(create_label(r'RMSD$\hspace{1}\theta$', ''))
    ax.set_title('(a)', loc='center')
    ax.set_ylim([0, 4])
    ax.legend()

    # (b)
    # sand and peat MB for radius
    ax = fig.add_subplot(4, 2, 2)
    ax.plot(rvec_sand, theta_nom_mb_vec_sand_radius, label='Heating and Cooling DP Sand')
    ax.plot(rvec_sand, theta_sig_mb_vec_sand_radius, label='Signal Processing SP and DP Sand')
    ax.plot(rvec_peat, theta_nom_mb_vec_peat_radius, label='Heating and Cooling DP Peat')
    ax.plot(rvec_peat, theta_sig_mb_vec_peat_radius, label='Signal Processing SP and DP Peat')
    ax.set_xlabel(create_label('$r$', 'mm'))
    ax.set_ylabel(create_label(r'MB$\hspace{1}\theta$', ''))
    ax.set_title('(b)', loc='center')

    # (c)
    # sand and peat RMSD for theta_o
    ax = fig.add_subplot(4, 2, 3)
    ax.plot(dvec_sand_theta_o, theta_nom_rmsd_vec_sand_theta_o)
    ax.plot(dvec_sand_theta_o, theta_sig_rmsd_vec_sand_theta_o)
    ax.plot(dvec_peat_theta_o, theta_nom_rmsd_vec_peat_theta_o)
    ax.plot(dvec_peat_theta_o, theta_sig_rmsd_vec_peat_theta_o)
    ax.set_xlabel(create_label(r'$\theta_o$', ''))
    ax.set_ylabel(create_label(r'RMSD$\hspace{1}\theta$', ''))
    ax.set_title('(c)', loc='center')

    # (d)
    # sand and peat MB for theta_o
    ax = fig.add_subplot(4, 2, 4)
    ax.plot(dvec_sand_theta_o, theta_nom_mb_vec_sand_theta_o)
    ax.plot(dvec_sand_theta_o, theta_sig_mb_vec_sand_theta_o)
    ax.plot(dvec_peat_theta_o, theta_nom_mb_vec_peat_theta_o)
    ax.plot(dvec_peat_theta_o, theta_sig_mb_vec_peat_theta_o)
    ax.set_xlabel(create_label(r'$\theta_o$', ''))
    ax.set_ylabel(create_label(r'MB$\hspace{1}\theta$', ''))
    ax.set_title('(d)', loc='center')

    # (e)
    # sand and peat RMSD for theta_m
    ax = fig.add_subplot(4, 2, 5)
    ax.plot(dvec_sand_theta_m, theta_nom_rmsd_vec_sand_theta_m)
    ax.plot(dvec_sand_theta_m, theta_sig_rmsd_vec_sand_theta_m)
    ax.plot(dvec_peat_theta_m, theta_nom_rmsd_vec_peat_theta_m)
    ax.plot(dvec_peat_theta_m, theta_sig_rmsd_vec_peat_theta_m)
    ax.set_xlabel(create_label(r'$\theta_m$', ''))
    ax.set_ylabel(create_label(r'RMSD$\hspace{1}\theta$', ''))
    ax.set_title('(e)', loc='center')

    # (f)
    # sand and peat MB for theta_m
    ax = fig.add_subplot(4, 2, 6)
    ax.plot(dvec_sand_theta_m, theta_nom_mb_vec_sand_theta_m)
    ax.plot(dvec_sand_theta_m, theta_sig_mb_vec_sand_theta_m)
    ax.plot(dvec_peat_theta_m, theta_nom_mb_vec_peat_theta_m)
    ax.plot(dvec_peat_theta_m, theta_sig_mb_vec_peat_theta_m)
    ax.set_xlabel(create_label(r'$\theta_m$', ''))
    ax.set_ylabel(create_label(r'MB$\hspace{1}\theta$', ''))
    ax.set_title('(f)', loc='center')

    # (g)
    # time
    ax = fig.add_subplot(4, 2, 7)
    # ax.plot(dvec_sand_time, theta_nom_rmsd_vec_sand_time)
    ax.plot(dvec_sand_time, theta_sig_rmsd_vec_sand_time, color=colors[1])
    # ax.plot(dvec_sand_time, theta_nom_rmsd_vec_peat_time)
    ax.plot(dvec_sand_time, theta_sig_rmsd_vec_peat_time, color=colors[3])
    ax.set_xlabel(create_label(r'$t_a$', 's'))
    ax.set_ylabel(create_label(r'RMSD$\hspace{1}\theta$', ''))
    ax.set_title('(g)', loc='center')

    # (h)
    ax = fig.add_subplot(4, 2, 8)
    # ax.plot(dvec_sand_time, theta_nom_mb_vec_sand_time)
    ax.plot(dvec_sand_time, theta_sig_mb_vec_sand_time, color=colors[1])
    # ax.plot(dvec_sand_time, theta_nom_mb_vec_peat_time)
    ax.plot(dvec_sand_time, theta_sig_mb_vec_peat_time, color=colors[3])
    ax.set_xlabel(create_label(r'$t_a$', 's'))
    ax.set_ylabel(create_label(r'MB$\hspace{1}\theta$', ''))
    ax.set_title('(h)', loc='center')

    plt.tight_layout()

    block = False
    plt.savefig(OAT_SENSITIVITY_FIG_NAME_THETA)
    if show:
        plt.show(block)

    #########################
    # RHO ONLY
    #########################

    fig = plt.figure(num=None, figsize=VLONG_FIGSIZE)

    # (a)
    ax = fig.add_subplot(4, 2, 1)
    ax.plot(rvec_sand, rho_nom_rmse_vec_sand_radius, label='Heating and Cooling DP Sand')
    ax.plot(rvec_sand, rho_sig_rmse_vec_sand_radius, label='Signal Processing SP and DP Sand')
    ax.plot(rvec_peat, rho_nom_rmse_vec_peat_radius, label='Heating and Cooling DP Peat')
    ax.plot(rvec_peat, rho_sig_rmse_vec_peat_radius, label='Signal Processing SP and DP Peat')
    ax.set_xlabel(create_label('$r$', 'mm'))
    ax.set_ylabel(create_label(r'RMSD$\hspace{1}\rho$', 'kg m^-3'))
    ax.set_title('(a)', loc='center')
    ax.legend()

    # (b)
    # sand and peat MB for radius
    ax = fig.add_subplot(4, 2, 2)
    ax.plot(rvec_sand, rho_nom_mb_vec_sand_radius, label='Heating and Cooling DP Sand')
    ax.plot(rvec_sand, rho_sig_mb_vec_sand_radius, label='Signal Processing SP and DP Sand')
    ax.plot(rvec_peat, rho_nom_mb_vec_peat_radius, label='Heating and Cooling DP Peat')
    ax.plot(rvec_peat, rho_sig_mb_vec_peat_radius, label='Signal Processing SP and DP Peat')
    ax.set_xlabel(create_label('$r$', 'mm'))
    ax.set_ylabel(create_label(r'MB$\hspace{1}\rho$', 'kg m^-3'))
    ax.set_title('(b)', loc='center')
    ax.set_ylim([None, 20000])

    # (c)
    # sand and peat RMSD for theta_o
    ax = fig.add_subplot(4, 2, 3)
    ax.plot(dvec_sand_theta_o, rho_nom_rmse_vec_sand_theta_o)
    ax.plot(dvec_sand_theta_o, rho_sig_rmse_vec_sand_theta_o)
    ax.plot(dvec_peat_theta_o, rho_nom_rmse_vec_peat_theta_o)
    ax.plot(dvec_peat_theta_o, rho_sig_rmse_vec_peat_theta_o)
    ax.set_xlabel(create_label(r'$\theta_o$', ''))
    ax.set_ylabel(create_label(r'RMSD$\hspace{1}\rho$', 'kg m^-3'))
    ax.set_title('(c)', loc='center')

    # (d)
    # sand and peat MB for theta_o
    ax = fig.add_subplot(4, 2, 4)
    ax.plot(dvec_sand_theta_o, rho_nom_mb_vec_sand_theta_o)
    ax.plot(dvec_sand_theta_o, rho_sig_mb_vec_sand_theta_o)
    ax.plot(dvec_peat_theta_o, rho_nom_mb_vec_peat_theta_o)
    ax.plot(dvec_peat_theta_o, rho_sig_mb_vec_peat_theta_o)
    ax.set_xlabel(create_label(r'$\theta_o$', ''))
    ax.set_ylabel(create_label(r'MB$\hspace{1}\rho$', 'kg m^-3'))
    ax.set_title('(d)', loc='center')

    # (e)
    # sand and peat RMSD for theta_m
    ax = fig.add_subplot(4, 2, 5)
    ax.plot(dvec_sand_theta_m, rho_nom_rmse_vec_sand_theta_m)
    ax.plot(dvec_sand_theta_m, rho_sig_rmse_vec_sand_theta_m)
    ax.plot(dvec_peat_theta_m, rho_nom_rmse_vec_peat_theta_m)
    ax.plot(dvec_peat_theta_m, rho_sig_rmse_vec_peat_theta_m)
    ax.set_xlabel(create_label(r'$\theta_m$', ''))
    ax.set_ylabel(create_label(r'RMSD$\hspace{1}\rho$', 'kg m^-3'))
    ax.set_title('(e)', loc='center')

    # (f)
    # sand and peat MB for theta_m
    ax = fig.add_subplot(4, 2, 6)
    ax.plot(dvec_sand_theta_m, rho_nom_mb_vec_sand_theta_m)
    ax.plot(dvec_sand_theta_m, rho_sig_mb_vec_sand_theta_m)
    ax.plot(dvec_peat_theta_m, rho_nom_mb_vec_peat_theta_m)
    ax.plot(dvec_peat_theta_m, rho_sig_mb_vec_peat_theta_m)
    ax.set_xlabel(create_label(r'$\theta_m$', ''))
    ax.set_ylabel(create_label(r'MB$\hspace{1}\rho$', 'kg m^-3'))
    ax.set_title('(f)', loc='center')

    # (g)
    # time
    ax = fig.add_subplot(4, 2, 7)
    # ax.plot(dvec_sand_time, rho_nom_rmse_vec_sand_time)
    ax.plot(dvec_sand_time, rho_sig_rmse_vec_sand_time, color=colors[1])
    # ax.plot(dvec_sand_time, rho_nom_rmse_vec_peat_time)
    ax.plot(dvec_sand_time, rho_sig_rmse_vec_peat_time, color=colors[3])
    ax.set_xlabel(create_label(r'$t_a$', 's'))
    ax.set_ylabel(create_label(r'RMSD$\hspace{1}\rho$', 'kg m^-3'))
    ax.set_title('(g)', loc='center')

    # (h)
    # time
    ax = fig.add_subplot(4, 2, 8)
    # ax.plot(dvec_sand_time, rho_nom_mb_vec_sand_time)
    ax.plot(dvec_sand_time, rho_sig_mb_vec_sand_time, color=colors[1])
    # ax.plot(dvec_sand_time, rho_nom_mb_vec_peat_time)
    ax.plot(dvec_sand_time, rho_sig_mb_vec_peat_time, color=colors[3])
    ax.set_xlabel(create_label(r'$t_a$', 's'))
    ax.set_ylabel(create_label(r'MB$\hspace{1}\rho$', 'kg m^-3'))
    ax.set_title('(h)', loc='center')

    plt.tight_layout()

    block = True
    plt.savefig(OAT_SENSITIVITY_FIG_NAME_RHO)
    if show:
        plt.show()
# DONE


def main():
    check = True                    # set check = True to check and see if the files exist (and do not re-run)
    # WARNING: The sensitivity analysis can take a number of hours
    run_all_process(check)          # run all processing to obtain the OAT files
    compute_rmsd_mb_oat()           # compute RMSD and MB for all of the data
    plot_oat(show=False)            # create the plots


if __name__ == '__main__':
    main()

