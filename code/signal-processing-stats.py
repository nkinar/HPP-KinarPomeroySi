import pickle
from get_size import length
from comparisons import compute_rmse, compute_mb, compute_percentage_diff, compute_variance, compute_sd, compute_pd
from constants import *
from PrettyTableWrapper import PrettyTableWrapper
from bs4 import BeautifulSoup
import pandas as pd


def cut_seq_trial(seq, start_num, end_num):
    """
    Cut the sequence belonging to a given trial
    :param seq:
    :param start_num:
    :param end_num:
    :return:
    """
    return seq[start_num-1:end_num]


def data_analysis_for_trials():
    """
    Perform the data analysis for all trials
    :return:
    """
    print('Loading data...')
    experiments = pickle.load(open(DATA_PICKLE_FILE_VEC, "rb"))

    theta_nom_sand, rho_nom_sand, theta_sig_sand, rho_sig_sand, \
    alpha_det_sand, alpha_sig_sand, k_det_sand, kdet_sig_sand, kdet_linear_sand, qav_vec_sand, qknown_vec_sand\
        = experiments[SAND_NAME]

    theta_nom_peat, rho_nom_peat, theta_sig_peat, rho_sig_peat, \
    alpha_det_peat, alpha_sig_peat, k_det_peat, kdet_sig_peat, kdet_linear_peat, qav_vec_peat, qknown_vec_peat\
        = experiments[PEAT_NAME]
    print('Done loading data...')

    ################################################################

    # compute the RMSD, MB and PD for sand and peat q (the heat input into the soil)
    print('Computing q comparisons')
    rmsd_sand_q = compute_rmse(qknown_vec_sand, qav_vec_sand)
    rmsd_peat_q = compute_rmse(qknown_vec_peat, qav_vec_peat)
    mb_sand_q = compute_mb(qknown_vec_sand, qav_vec_sand)
    mb_peat_q = compute_mb(qknown_vec_peat, qav_vec_peat)
    pd_sand_q = compute_pd(qknown_vec_sand, qav_vec_sand)
    pd_peat_q = compute_pd(qknown_vec_peat, qav_vec_peat)
    comp_list = [('sand',  rmsd_sand_q, mb_sand_q,  pd_sand_q),
                 ('peat', rmsd_peat_q, mb_peat_q, pd_peat_q)]
    labels = ['soil', 'rmsd', 'mb', 'pd']
    df_q = pd.DataFrame.from_records(comp_list, columns=labels)
    df_q.to_html(TABLE_PATH + TABLE_NAME_Q_COMP + HTML_EXT)
    print('DONE saving out HTML comparisons for q')

    ################################################################

    # first sequence of experiments
    ns = length(theta_nom_sand)
    xs = np.linspace(1, ns, ns)

    # second sequence of experiments
    nnp = length(theta_nom_peat)
    xp = np.linspace(xs[-1]+1, xs[-1]+nnp, nnp)

    # break the sequences up for each set so that the stats can be run

    ########################
    # NUMBERS SAND
    ########################
    start_num = 1
    end_num = 5
    calc_sand_theta_nom_sand_num_1_5 = cut_seq_trial(theta_nom_sand, start_num, end_num)
    calc_sand_rho_nom_sand_num_1_5 = cut_seq_trial(rho_nom_sand, start_num, end_num)
    calc_sand_theta_sig_sand_num_1_5 = cut_seq_trial(theta_sig_sand, start_num, end_num)
    calc_sand_rho_sig_sand_num_1_5 = cut_seq_trial(rho_sig_sand, start_num, end_num)

    # theta nominal for sand 1-5
    rmse__calc_sand_theta_nom_sand_num_1_5 = compute_rmse(THETA_NOMINAL_SAND, calc_sand_theta_nom_sand_num_1_5)
    mb__calc_sand_theta_nom_sand_num_1_5 = compute_mb(THETA_NOMINAL_SAND, calc_sand_theta_nom_sand_num_1_5)
    pd__calc_sand_theta_nom_sand_num_1_5 = compute_percentage_diff(THETA_NOMINAL_SAND, calc_sand_theta_nom_sand_num_1_5)
    variance__calc_sand_theta_nom_sand_num_1_5 = compute_variance(calc_sand_theta_nom_sand_num_1_5)
    sd__calc_sand_theta_nom_sand_num_1_5 = compute_sd(calc_sand_theta_nom_sand_num_1_5)

    # rho nominal for sand 1-5
    rmse__calc_sand_rho_nom_sand_num_1_5 = compute_rmse(DENSITY_NOMINAL_SAND, calc_sand_rho_nom_sand_num_1_5)
    mb__calc_sand_rho_nom_sand_num_1_5 = compute_mb(DENSITY_NOMINAL_SAND, calc_sand_rho_nom_sand_num_1_5)
    pd__calc_sand_rho_nom_sand_num_1_5 = compute_percentage_diff(DENSITY_NOMINAL_SAND, calc_sand_rho_nom_sand_num_1_5)
    variance__calc_rho_theta_nom_sand_num_1_5 = compute_variance(calc_sand_rho_nom_sand_num_1_5)
    sd__calc_sand_rho_nom_sand_num_1_5 = compute_sd(calc_sand_rho_nom_sand_num_1_5)

    # theta signal for sand 1-5
    rmse__calc_sand_theta_sig_sand_num_1_5 = compute_rmse(THETA_NOMINAL_SAND, calc_sand_theta_sig_sand_num_1_5)
    mb__calc_sand_theta_sig_sand_num_1_5 = compute_mb(THETA_NOMINAL_SAND, calc_sand_theta_sig_sand_num_1_5)
    pd__calc_sand_theta_sig_sand_num_1_5 = compute_percentage_diff(THETA_NOMINAL_SAND, calc_sand_theta_sig_sand_num_1_5)
    variance__calc_sand_theta_sig_sand_num_1_5 = compute_variance(calc_sand_theta_sig_sand_num_1_5)
    sd__calc_sand_theta_sig_sand_num_1_5 = compute_sd(calc_sand_theta_sig_sand_num_1_5)

    # rho signal for sand 1-5
    rmse__calc_sand_rho_sig_sand_num_1_5 = compute_rmse(DENSITY_NOMINAL_SAND, calc_sand_rho_sig_sand_num_1_5)
    mb__calc_sand_rho_sig_sand_num_1_5 = compute_mb(DENSITY_NOMINAL_SAND, calc_sand_rho_sig_sand_num_1_5)
    pd__calc_sand_rho_sig_sand_num_1_5 = compute_percentage_diff(DENSITY_NOMINAL_SAND, calc_sand_rho_sig_sand_num_1_5)
    variance__calc_rho_sig_nom_sand_num_1_5 = compute_variance(calc_sand_rho_sig_sand_num_1_5)
    sd__calc_sand_rho_sig_sand_num_1_5 = compute_sd(calc_sand_rho_sig_sand_num_1_5)

    ##########################################################################################

    start_num = 6
    end_num = 10
    calc_sand_theta_nom_sand_num_6_10 = cut_seq_trial(theta_nom_sand, start_num, end_num)
    calc_sand_rho_nom_sand_num_6_10 = cut_seq_trial(rho_nom_sand, start_num, end_num)
    calc_sand_theta_sig_sand_num_6_10 = cut_seq_trial(theta_sig_sand, start_num, end_num)
    calc_sand_rho_sig_sand_num_6_10 = cut_seq_trial(rho_sig_sand, start_num, end_num)

    # theta nominal for sand 6-10
    rmse__calc_sand_theta_nom_sand_num_6_10 = compute_rmse(THETA_NOMINAL_SAND, calc_sand_theta_nom_sand_num_6_10)
    mb__calc_sand_theta_nom_sand_num_6_10 = compute_mb(THETA_NOMINAL_SAND, calc_sand_theta_nom_sand_num_6_10)
    pd__calc_sand_theta_nom_sand_num_6_10 = compute_percentage_diff(THETA_NOMINAL_SAND, calc_sand_theta_nom_sand_num_6_10)
    variance__calc_sand_theta_nom_sand_num_6_10 = compute_variance(calc_sand_theta_nom_sand_num_6_10)
    sd__calc_sand_theta_nom_sand_num_6_10 = compute_sd(calc_sand_theta_nom_sand_num_6_10)

    # rho nominal for sand 6-10
    rmse__calc_sand_rho_nom_sand_num_6_10 = compute_rmse(DENSITY_NOMINAL_SAND, calc_sand_rho_nom_sand_num_6_10)
    mb__calc_sand_rho_nom_sand_num_6_10 = compute_mb(DENSITY_NOMINAL_SAND, calc_sand_rho_nom_sand_num_6_10)
    pd__calc_sand_rho_nom_sand_num_6_10 = compute_percentage_diff(DENSITY_NOMINAL_SAND, calc_sand_rho_nom_sand_num_6_10)
    variance__calc_rho_theta_nom_sand_num_6_10 = compute_variance(calc_sand_rho_nom_sand_num_6_10)
    sd__calc_sand_rho_nom_sand_num_6_10 = compute_sd(calc_sand_rho_nom_sand_num_6_10)

    # theta signal for sand 6-10
    rmse__calc_sand_theta_sig_sand_num_6_10 = compute_rmse(THETA_NOMINAL_SAND, calc_sand_theta_sig_sand_num_6_10)
    mb__calc_sand_theta_sig_sand_num_6_10 = compute_mb(THETA_NOMINAL_SAND, calc_sand_theta_sig_sand_num_6_10)
    pd__calc_sand_theta_sig_sand_num_6_10 = compute_percentage_diff(THETA_NOMINAL_SAND, calc_sand_theta_sig_sand_num_6_10)
    variance__calc_sand_theta_sig_sand_num_6_10 = compute_variance(calc_sand_theta_sig_sand_num_6_10)
    sd__calc_sand_theta_sig_sand_num_6_10 = compute_sd(calc_sand_theta_sig_sand_num_6_10)

    # rho signal for sand 6-10
    rmse__calc_sand_rho_sig_sand_num_6_10 = compute_rmse(DENSITY_NOMINAL_SAND, calc_sand_rho_sig_sand_num_6_10)
    mb__calc_sand_rho_sig_sand_num_6_10 = compute_mb(DENSITY_NOMINAL_SAND, calc_sand_rho_sig_sand_num_6_10)
    pd__calc_sand_rho_sig_sand_num_6_10 = compute_percentage_diff(DENSITY_NOMINAL_SAND, calc_sand_rho_sig_sand_num_6_10)
    variance__calc_rho_sig_nom_sand_num_6_10 = compute_variance(calc_sand_rho_sig_sand_num_6_10)
    sd__calc_sand_rho_sig_sand_num_6_10 = compute_sd(calc_sand_rho_sig_sand_num_6_10)

    ##########################################################################################

    start_num = 11
    end_num = 15
    calc_sand_theta_nom_sand_num_11_15 = cut_seq_trial(theta_nom_sand, start_num, end_num)
    calc_sand_rho_nom_sand_num_11_15 = cut_seq_trial(rho_nom_sand, start_num, end_num)
    calc_sand_theta_sig_sand_num_11_15 = cut_seq_trial(theta_sig_sand, start_num, end_num)
    calc_sand_rho_sig_sand_num_11_15 = cut_seq_trial(rho_sig_sand, start_num, end_num)

    # theta nominal for sand 11-15
    rmse__calc_sand_theta_nom_sand_num_11_15 = compute_rmse(THETA_NOMINAL_SAND, calc_sand_theta_nom_sand_num_11_15)
    mb__calc_sand_theta_nom_sand_num_11_15 = compute_mb(THETA_NOMINAL_SAND, calc_sand_theta_nom_sand_num_11_15)
    pd__calc_sand_theta_nom_sand_num_11_15 = compute_percentage_diff(THETA_NOMINAL_SAND, calc_sand_theta_nom_sand_num_11_15)
    variance__calc_sand_theta_nom_sand_num_11_15 = compute_variance(calc_sand_theta_nom_sand_num_11_15)
    sd__calc_sand_theta_nom_sand_num_11_15 = compute_sd(calc_sand_theta_nom_sand_num_11_15)

    # rho nominal for sand 11-15
    rmse__calc_sand_rho_nom_sand_num_11_15 = compute_rmse(DENSITY_NOMINAL_SAND, calc_sand_rho_nom_sand_num_11_15)
    mb__calc_sand_rho_nom_sand_num_11_15 = compute_mb(DENSITY_NOMINAL_SAND, calc_sand_rho_nom_sand_num_11_15)
    pd__calc_sand_rho_nom_sand_num_11_15 = compute_percentage_diff(DENSITY_NOMINAL_SAND, calc_sand_rho_nom_sand_num_11_15)
    variance__calc_rho_theta_nom_sand_num_11_15 = compute_variance(calc_sand_rho_nom_sand_num_11_15)
    sd__calc_sand_rho_nom_sand_num_11_15 = compute_sd(calc_sand_rho_nom_sand_num_11_15)

    # theta signal for sand 11-15
    rmse__calc_sand_theta_sig_sand_num_11_15 = compute_rmse(THETA_NOMINAL_SAND, calc_sand_theta_sig_sand_num_11_15)
    mb__calc_sand_theta_sig_sand_num_11_15 = compute_mb(THETA_NOMINAL_SAND, calc_sand_theta_sig_sand_num_11_15)
    pd__calc_sand_theta_sig_sand_num_11_15 = compute_percentage_diff(THETA_NOMINAL_SAND, calc_sand_theta_sig_sand_num_11_15)
    variance__calc_sand_theta_sig_sand_num_11_15 = compute_variance(calc_sand_theta_sig_sand_num_11_15)
    sd__calc_sand_theta_sig_sand_num_11_15 = compute_sd(calc_sand_theta_sig_sand_num_11_15)

    # rho signal for sand 11-15
    rmse__calc_sand_rho_sig_sand_num_11_15 = compute_rmse(DENSITY_NOMINAL_SAND, calc_sand_rho_sig_sand_num_11_15)
    mb__calc_sand_rho_sig_sand_num_11_15 = compute_mb(DENSITY_NOMINAL_SAND, calc_sand_rho_sig_sand_num_11_15)
    pd__calc_sand_rho_sig_sand_num_11_15 = compute_percentage_diff(DENSITY_NOMINAL_SAND, calc_sand_rho_sig_sand_num_11_15)
    variance__calc_rho_sig_nom_sand_num_11_15 = compute_variance(calc_sand_rho_sig_sand_num_11_15)
    sd__calc_sand_rho_sig_sand_num_11_15 = compute_sd(calc_sand_rho_sig_sand_num_11_15)

    #################################################################################################

    start_num = 16
    end_num = 20
    calc_sand_theta_nom_sand_num_16_20 = cut_seq_trial(theta_nom_sand, start_num, end_num)
    calc_sand_rho_nom_sand_num_16_20 = cut_seq_trial(rho_nom_sand, start_num, end_num)
    calc_sand_theta_sig_sand_num_16_20 = cut_seq_trial(theta_sig_sand, start_num, end_num)
    calc_sand_rho_sig_sand_num_16_20 = cut_seq_trial(rho_sig_sand, start_num, end_num)

    # theta nominal for sand 16-20
    rmse__calc_sand_theta_nom_sand_num_16_20 = compute_rmse(THETA_NOMINAL_SAND, calc_sand_theta_nom_sand_num_16_20)
    mb__calc_sand_theta_nom_sand_num_16_20 = compute_mb(THETA_NOMINAL_SAND, calc_sand_theta_nom_sand_num_16_20)
    pd__calc_sand_theta_nom_sand_num_16_20 = compute_percentage_diff(THETA_NOMINAL_SAND,
                                                                     calc_sand_theta_nom_sand_num_16_20)
    variance__calc_sand_theta_nom_sand_num_16_20 = compute_variance(calc_sand_theta_nom_sand_num_16_20)
    sd__calc_sand_theta_nom_sand_num_16_20 = compute_sd(calc_sand_theta_nom_sand_num_16_20)

    # rho nominal for sand 16-20
    rmse__calc_sand_rho_nom_sand_num_16_20 = compute_rmse(DENSITY_NOMINAL_SAND, calc_sand_rho_nom_sand_num_16_20)
    mb__calc_sand_rho_nom_sand_num_16_20 = compute_mb(DENSITY_NOMINAL_SAND, calc_sand_rho_nom_sand_num_16_20)
    pd__calc_sand_rho_nom_sand_num_16_20 = compute_percentage_diff(DENSITY_NOMINAL_SAND,
                                                                   calc_sand_rho_nom_sand_num_16_20)
    variance__calc_rho_theta_nom_sand_num_16_20 = compute_variance(calc_sand_rho_nom_sand_num_16_20)
    sd__calc_sand_rho_nom_sand_num_16_20 = compute_sd(calc_sand_rho_nom_sand_num_16_20)

    # theta signal for sand 16-20
    rmse__calc_sand_theta_sig_sand_num_16_20 = compute_rmse(THETA_NOMINAL_SAND, calc_sand_theta_sig_sand_num_16_20)
    mb__calc_sand_theta_sig_sand_num_16_20 = compute_mb(THETA_NOMINAL_SAND, calc_sand_theta_sig_sand_num_16_20)
    pd__calc_sand_theta_sig_sand_num_16_20 = compute_percentage_diff(THETA_NOMINAL_SAND,
                                                                     calc_sand_theta_sig_sand_num_16_20)
    variance__calc_sand_theta_sig_sand_num_16_20 = compute_variance(calc_sand_theta_sig_sand_num_16_20)
    sd__calc_sand_theta_sig_sand_num_16_20 = compute_sd(calc_sand_theta_sig_sand_num_16_20)

    # rho signal for sand 16-20
    rmse__calc_sand_rho_sig_sand_num_16_20 = compute_rmse(DENSITY_NOMINAL_SAND, calc_sand_rho_sig_sand_num_16_20)
    mb__calc_sand_rho_sig_sand_num_16_20 = compute_mb(DENSITY_NOMINAL_SAND, calc_sand_rho_sig_sand_num_16_20)
    pd__calc_sand_rho_sig_sand_num_16_20 = compute_percentage_diff(DENSITY_NOMINAL_SAND,
                                                                   calc_sand_rho_sig_sand_num_16_20)
    variance__calc_rho_sig_nom_sand_num_16_20 = compute_variance(calc_sand_rho_sig_sand_num_16_20)
    sd__calc_sand_rho_sig_sand_num_16_20 = compute_sd(calc_sand_rho_sig_sand_num_16_20)

    ########################
    # NUMBERS PEAT
    ########################

    start_num = 1
    end_num = 5
    calc_peat_theta_nom_peat_num_21_25 = cut_seq_trial(theta_nom_peat, start_num, end_num)
    calc_peat_rho_nom_peat_num_21_25 = cut_seq_trial(rho_nom_peat, start_num, end_num)
    calc_peat_theta_sig_peat_num_21_25 = cut_seq_trial(theta_sig_peat, start_num, end_num)
    calc_peat_rho_sig_peat_num_21_25 = cut_seq_trial(rho_sig_peat, start_num, end_num)

    # theta nominal for peat 21-25
    rmse__calc_peat_theta_nom_peat_num_21_25 = compute_rmse(THETA_NOMINAL_PEAT, calc_peat_theta_nom_peat_num_21_25)
    mb__calc_peat_theta_nom_peat_num_21_25 = compute_mb(THETA_NOMINAL_PEAT, calc_peat_theta_nom_peat_num_21_25)
    pd__calc_peat_theta_nom_peat_num_21_25 = compute_percentage_diff(THETA_NOMINAL_PEAT, calc_peat_theta_nom_peat_num_21_25)
    variance__calc_peat_theta_nom_peat_num_21_25 = compute_variance(calc_peat_theta_nom_peat_num_21_25)
    sd__calc_peat_theta_nom_peat_num_21_25 = compute_sd(calc_peat_theta_nom_peat_num_21_25)

    # rho nominal for peat 21-25
    rmse__calc_peat_rho_nom_sand_num_21_25 = compute_rmse(DENSITY_NOMINAL_PEAT, calc_peat_rho_nom_peat_num_21_25)
    mb__calc_peat_rho_nom_sand_num_21_25 = compute_mb(DENSITY_NOMINAL_PEAT, calc_peat_rho_nom_peat_num_21_25)
    pd__calc_peat_rho_nom_sand_num_21_25 = compute_percentage_diff(DENSITY_NOMINAL_PEAT, calc_peat_rho_nom_peat_num_21_25)
    variance__calc_peat_rho_theta_nom_sand_num_21_25 = compute_variance(calc_peat_rho_nom_peat_num_21_25)
    sd__calc_peat_rho_nom_sand_num_21_25 = compute_sd(calc_peat_rho_nom_peat_num_21_25)

    # theta signal for peat 21-25
    rmse__calc_peat_theta_sig_sand_num_21_25 = compute_rmse(THETA_NOMINAL_PEAT, calc_peat_theta_sig_peat_num_21_25)
    mb__calc_peat_theta_sig_sand_num_21_25 = compute_mb(THETA_NOMINAL_PEAT, calc_peat_theta_sig_peat_num_21_25)
    pd__calc_peat_theta_sig_sand_num_21_25 = compute_percentage_diff(THETA_NOMINAL_PEAT, calc_peat_theta_sig_peat_num_21_25)
    variance__calc_peat_theta_sig_sand_num_21_25 = compute_variance(calc_peat_theta_sig_peat_num_21_25)
    sd__calc_peat_theta_sig_sand_num_21_25 = compute_sd(calc_peat_theta_sig_peat_num_21_25)

    # rho signal for peat 21-25
    rmse__calc_peat_rho_sig_peat_num_21_25 = compute_rmse(DENSITY_NOMINAL_PEAT, calc_peat_rho_sig_peat_num_21_25)
    mb__calc_peat_rho_sig_peat_num_21_25 = compute_mb(DENSITY_NOMINAL_PEAT, calc_peat_rho_sig_peat_num_21_25)
    pd__calc_peat_rho_sig_peat_num_21_25 = compute_percentage_diff(DENSITY_NOMINAL_PEAT, calc_peat_rho_sig_peat_num_21_25)
    variance__calc_peat_rho_sig_nom_peat_num_21_25 = compute_variance(calc_peat_rho_sig_peat_num_21_25)
    sd__calc_peat_rho_sig_peat_num_21_25 = compute_sd(calc_peat_rho_sig_peat_num_21_25)

    #########################################################################################################

    start_num = 6
    end_num = 10
    calc_peat_theta_nom_peat_num_26_30 = cut_seq_trial(theta_nom_peat, start_num, end_num)
    calc_peat_rho_nom_peat_num_26_30 = cut_seq_trial(rho_nom_peat, start_num, end_num)
    calc_peat_theta_sig_peat_num_26_30 = cut_seq_trial(theta_sig_peat, start_num, end_num)
    calc_peat_rho_sig_peat_num_26_30 = cut_seq_trial(rho_sig_peat, start_num, end_num)

    # theta nominal for peat 26-30
    rmse__calc_peat_theta_nom_peat_num_26_30 = compute_rmse(THETA_NOMINAL_PEAT, calc_peat_theta_nom_peat_num_26_30)
    mb__calc_peat_theta_nom_peat_num_26_30 = compute_mb(THETA_NOMINAL_PEAT, calc_peat_theta_nom_peat_num_26_30)
    pd__calc_peat_theta_nom_peat_num_26_30 = compute_percentage_diff(THETA_NOMINAL_PEAT, calc_peat_theta_nom_peat_num_26_30)
    variance__calc_peat_theta_nom_peat_num_26_30 = compute_variance(calc_peat_theta_nom_peat_num_26_30)
    sd__calc_peat_theta_nom_peat_num_26_30 = compute_sd(calc_peat_theta_nom_peat_num_26_30)

    # rho nominal for peat 26-30
    rmse__calc_peat_rho_nom_peat_num_26_30 = compute_rmse(DENSITY_NOMINAL_PEAT, calc_peat_rho_nom_peat_num_26_30)
    mb__calc_peat_rho_nom_peat_num_26_30 = compute_mb(DENSITY_NOMINAL_PEAT, calc_peat_rho_nom_peat_num_26_30)
    pd__calc_peat_rho_nom_peat_num_26_30 = compute_percentage_diff(DENSITY_NOMINAL_PEAT, calc_peat_rho_nom_peat_num_26_30)
    variance__calc_peat_rho_theta_nom_peat_num_26_30 = compute_variance(calc_peat_rho_nom_peat_num_26_30)
    sd__calc_peat_rho_nom_peat_num_26_30 = compute_sd(calc_peat_rho_nom_peat_num_26_30)

    # theta signal for peat 26-30
    rmse__calc_peat_theta_sig_peat_num_26_30 = compute_rmse(THETA_NOMINAL_PEAT, calc_peat_theta_sig_peat_num_26_30)
    mb__calc_peat_theta_sig_peat_num_26_30 = compute_mb(THETA_NOMINAL_PEAT, calc_peat_theta_sig_peat_num_26_30)
    pd__calc_peat_theta_sig_peat_num_26_30 = compute_percentage_diff(THETA_NOMINAL_PEAT, calc_peat_theta_sig_peat_num_26_30)
    variance__calc_peat_theta_sig_peat_num_26_30 = compute_variance(calc_peat_theta_sig_peat_num_26_30)
    sd__calc_peat_theta_sig_peat_num_26_30 = compute_sd(calc_peat_theta_sig_peat_num_26_30)

    # rho signal for peat 26-30
    rmse__calc_peat_rho_sig_peat_num_26_30 = compute_rmse(DENSITY_NOMINAL_PEAT, calc_peat_rho_sig_peat_num_26_30)
    mb__calc_peat_rho_sig_peat_num_26_30 = compute_mb(DENSITY_NOMINAL_PEAT, calc_peat_rho_sig_peat_num_26_30)
    pd__calc_peat_rho_sig_peat_num_26_30 = compute_percentage_diff(DENSITY_NOMINAL_PEAT, calc_peat_rho_sig_peat_num_26_30)
    variance__calc_peat_rho_sig_nom_peat_num_26_30 = compute_variance(calc_peat_rho_sig_peat_num_26_30)
    sd__calc_peat_rho_sig_peat_num_26_30 = compute_sd(calc_peat_rho_sig_peat_num_26_30)

    #########################################################################################################

    start_num = 11
    end_num = 15
    calc_peat_theta_nom_peat_num_31_35 = cut_seq_trial(theta_nom_peat, start_num, end_num)
    calc_peat_rho_nom_peat_num_31_35 = cut_seq_trial(rho_nom_peat, start_num, end_num)
    calc_peat_theta_sig_peat_num_31_35 = cut_seq_trial(theta_sig_peat, start_num, end_num)
    calc_peat_rho_sig_peat_num_31_35 = cut_seq_trial(rho_sig_peat, start_num, end_num)

    # theta nominal for peat 31-25
    rmse__calc_peat_theta_nom_peat_num_31_35 = compute_rmse(THETA_NOMINAL_PEAT, calc_peat_theta_nom_peat_num_31_35)
    mb__calc_peat_theta_nom_peat_num_31_35 = compute_mb(THETA_NOMINAL_PEAT, calc_peat_theta_nom_peat_num_31_35)
    pd__calc_peat_theta_nom_peat_num_31_35 = compute_percentage_diff(THETA_NOMINAL_PEAT, calc_peat_theta_nom_peat_num_31_35)
    variance__calc_peat_theta_nom_peat_num_31_35 = compute_variance(calc_peat_theta_nom_peat_num_31_35)
    sd__calc_peat_theta_nom_peat_num_31_35 = compute_sd(calc_peat_theta_nom_peat_num_31_35)

    # rho nominal for peat 31-25
    rmse__calc_peat_rho_nom_peat_num_31_35 = compute_rmse(DENSITY_NOMINAL_PEAT, calc_peat_rho_nom_peat_num_31_35)
    mb__calc_peat_rho_nom_peat_num_31_35 = compute_mb(DENSITY_NOMINAL_PEAT, calc_peat_rho_nom_peat_num_31_35)
    pd__calc_peat_rho_nom_peat_num_31_35 = compute_percentage_diff(DENSITY_NOMINAL_PEAT, calc_peat_rho_nom_peat_num_31_35)
    variance__calc_peat_rho_theta_nom_peat_num_31_35 = compute_variance(calc_peat_rho_nom_peat_num_31_35)
    sd__calc_peat_rho_nom_peat_num_31_35 = compute_sd(calc_peat_rho_nom_peat_num_31_35)

    # theta signal for peat 31-25
    rmse__calc_peat_theta_sig_peat_num_31_35 = compute_rmse(THETA_NOMINAL_PEAT, calc_peat_theta_sig_peat_num_31_35)
    mb__calc_peat_theta_sig_peat_num_31_35 = compute_mb(THETA_NOMINAL_PEAT, calc_peat_theta_sig_peat_num_31_35)
    pd__calc_peat_theta_sig_peat_num_31_35 = compute_percentage_diff(THETA_NOMINAL_PEAT, calc_peat_theta_sig_peat_num_31_35)
    variance__calc_peat_theta_sig_peat_num_31_35 = compute_variance(calc_peat_theta_sig_peat_num_31_35)
    sd__calc_peat_theta_sig_peat_num_31_35 = compute_sd(calc_peat_theta_sig_peat_num_31_35)

    # rho signal for peat 31-25
    rmse__calc_peat_rho_sig_peat_num_31_35 = compute_rmse(DENSITY_NOMINAL_PEAT, calc_peat_rho_sig_peat_num_31_35)
    mb__calc_peat_rho_sig_peat_num_31_35 = compute_mb(DENSITY_NOMINAL_PEAT, calc_peat_rho_sig_peat_num_31_35)
    pd__calc_peat_rho_sig_peat_num_31_35 = compute_percentage_diff(DENSITY_NOMINAL_PEAT, calc_peat_rho_sig_peat_num_31_35)
    variance__calc__peat_rho_sig_nom_peat_num_31_35 = compute_variance(calc_peat_rho_sig_peat_num_31_35)
    sd__calc_peat_rho_sig_peat_num_31_35 = compute_sd(calc_peat_rho_sig_peat_num_31_35)

    #########################################################################################################

    ###########################
    # ALL SAND
    ###########################

    calc_sand_theta_nom_sand_all = theta_nom_sand
    calc_sand_rho_nom_sand_all = rho_nom_sand
    calc_sand_theta_sig_sand_all = theta_sig_sand
    calc_sand_rho_sig_sand_all = rho_sig_sand

    rmse__calc_sand_theta_nom_sand_all = compute_rmse(THETA_NOMINAL_SAND, calc_sand_theta_nom_sand_all)
    mb__calc_sand_theta_nom_sand_all = compute_mb(THETA_NOMINAL_SAND, calc_sand_theta_nom_sand_all)
    pd__calc_sand_theta_nom_sand_all = compute_percentage_diff(THETA_NOMINAL_SAND, calc_sand_theta_nom_sand_all)
    variance__calc_sand_theta_nom_sand_all = compute_variance(calc_sand_theta_nom_sand_all)
    sd__calc_sand_theta_nom_sand_trial_all = compute_sd(calc_sand_theta_nom_sand_all)

    rmse__calc_sand_rho_nom_sand_all = compute_rmse(DENSITY_NOMINAL_SAND, calc_sand_rho_nom_sand_all)
    mb__calc_sand_rho_nom_sand_all = compute_mb(DENSITY_NOMINAL_SAND, calc_sand_rho_nom_sand_all)
    pd__calc_sand_rho_nom_sand_all = compute_percentage_diff(DENSITY_NOMINAL_SAND, calc_sand_rho_nom_sand_all)
    variance__calc_sand_rho_nom_sand_all = compute_variance(calc_sand_rho_nom_sand_all)
    sd__calc_sand_rho_nom_sand_all = compute_sd(calc_sand_rho_nom_sand_all)

    rmse__calc_sand_theta_sig_sand_all = compute_rmse(THETA_NOMINAL_SAND, calc_sand_theta_sig_sand_all)
    mb__calc_sand_theta_sig_sand_all = compute_mb(THETA_NOMINAL_SAND, calc_sand_theta_sig_sand_all)
    pd__calc_sand_theta_sig_sand_all = compute_percentage_diff(THETA_NOMINAL_SAND, calc_sand_theta_sig_sand_all)
    variance__calc_sand_theta_sig_sand_all = compute_variance(calc_sand_theta_sig_sand_all)
    sd__calc_sand_theta_sig_sand_all = compute_sd(calc_sand_theta_sig_sand_all)

    rmse__calc_sand_rho_sig_sand_all = compute_rmse(DENSITY_NOMINAL_SAND, calc_sand_rho_sig_sand_all)
    mb__calc_sand_rho_sig_sand_all = compute_mb(DENSITY_NOMINAL_SAND, calc_sand_rho_sig_sand_all)
    pd__calc_sand_rho_sig_sand_all = compute_percentage_diff(DENSITY_NOMINAL_SAND, calc_sand_rho_sig_sand_all)
    variance__calc_sand_rho_sig_sand_all = compute_variance(calc_sand_rho_sig_sand_all)
    sd__calc_sand_rho_sig_sand_all = compute_sd(calc_sand_rho_sig_sand_all)

    ###########################
    # ALL PEAT
    ###########################
    calc_peat_theta_nom_peat_all = theta_nom_peat
    calc_peat_rho_nom_peat_all = rho_nom_peat
    calc_peat_theta_sig_peat_all = theta_sig_peat
    calc_peat_rho_sig_peat_all = rho_sig_peat

    rmse__calc_peat_theta_nom_peat_all = compute_rmse(THETA_NOMINAL_PEAT, calc_peat_theta_nom_peat_all)
    mb__calc_peat_theta_nom_peat_all = compute_mb(THETA_NOMINAL_PEAT, calc_peat_theta_nom_peat_all)
    pd__calc_peat_theta_nom_peat_all = compute_percentage_diff(THETA_NOMINAL_PEAT, calc_peat_theta_nom_peat_all)
    variance__calc_peat_theta_nom_peat_all = compute_variance(calc_peat_theta_nom_peat_all)
    sd__calc_peat_theta_nom_peat_trial_all = compute_sd(calc_peat_theta_nom_peat_all)

    rmse__calc_peat_rho_nom_peat_all = compute_rmse(DENSITY_NOMINAL_PEAT, calc_peat_rho_nom_peat_all)
    mb__calc_peat_rho_nom_peat_all = compute_mb(DENSITY_NOMINAL_PEAT, calc_peat_rho_nom_peat_all)
    pd__calc_peat_rho_nom_peat_all = compute_percentage_diff(DENSITY_NOMINAL_PEAT, calc_peat_rho_nom_peat_all)
    variance__calc_peat_rho_nom_peat_all = compute_variance(calc_peat_rho_nom_peat_all)
    sd__calc_peat_rho_nom_peat_all = compute_sd(calc_peat_rho_nom_peat_all)

    rmse__calc_peat_theta_sig_peat_all = compute_rmse(THETA_NOMINAL_PEAT, calc_peat_theta_sig_peat_all)
    mb__calc_peat_theta_sig_peat_all = compute_mb(THETA_NOMINAL_PEAT, calc_peat_theta_sig_peat_all)
    pd__calc_peat_theta_sig_peat_all = compute_percentage_diff(THETA_NOMINAL_PEAT, calc_peat_theta_sig_peat_all)
    variance__calc_peat_theta_sig_peat_all = compute_variance(calc_peat_theta_sig_peat_all)
    sd__calc_peat_theta_sig_peat_all = compute_sd(calc_peat_theta_sig_peat_all)

    rmse__calc_peat_rho_sig_peat_all = compute_rmse(DENSITY_NOMINAL_PEAT, calc_peat_rho_sig_peat_all)
    mb__calc_peat_rho_sig_peat_all = compute_mb(DENSITY_NOMINAL_PEAT, calc_peat_rho_sig_peat_all)
    pd__calc_peat_rho_sig_peat_all = compute_percentage_diff(DENSITY_NOMINAL_PEAT, calc_peat_rho_sig_peat_all)
    variance__calc_peat_rho_sig_peat_all = compute_variance(calc_peat_rho_sig_peat_all)
    sd__calc_peat_rho_sig_peat_all = compute_sd(calc_peat_rho_sig_peat_all)

    # EXPORT ALL OF THE ABOVE IN AN HTML TABLE
    tabd = PrettyTableWrapper(DEC_PLACES_TAB)
    # 22 columns (list vertically)
    tabd.field_names(["Soil Type",
                        "Identifier",

                        "theta_nominal RMSE",
                        "theta_nominal MB",
                        "theta_nominal PD (%)",
                        "theta_nominal Variance",
                        "theta_nominal SD",

                        "rho_nominal RMSE",
                        "rho_nominal MB",
                        "rho_nominal PD (%)",
                        "rho_nominal Variance",
                        "rho_nominal SD",

                        "theta_sig RMSE",
                        "theta_sig MB",
                        "theta_sig PD (%)",
                        "theta_sig Variance",
                        "theta_sig SD",

                        "rho_sig RMSE",
                        "rho_sig MB",
                        "rho_sig PD (%)",
                        "rho_sig Variance",
                        "rho_sig SD"
                        ])

    tabd.add_row(['sand',
                  '1-5',
                  rmse__calc_sand_theta_nom_sand_num_1_5,
                  mb__calc_sand_theta_nom_sand_num_1_5,
                  pd__calc_sand_theta_nom_sand_num_1_5,
                  variance__calc_sand_theta_nom_sand_num_1_5,
                  sd__calc_sand_theta_nom_sand_num_1_5,

                  rmse__calc_sand_rho_nom_sand_num_1_5,
                  mb__calc_sand_rho_nom_sand_num_1_5,
                  pd__calc_sand_rho_nom_sand_num_1_5,
                  variance__calc_rho_theta_nom_sand_num_1_5,
                  sd__calc_sand_rho_nom_sand_num_1_5,

                  rmse__calc_sand_theta_sig_sand_num_1_5,
                  mb__calc_sand_theta_sig_sand_num_1_5,
                  pd__calc_sand_theta_sig_sand_num_1_5,
                  variance__calc_sand_theta_sig_sand_num_1_5,
                  sd__calc_sand_theta_sig_sand_num_1_5,

                  rmse__calc_sand_rho_sig_sand_num_1_5,
                  mb__calc_sand_rho_sig_sand_num_1_5,
                  pd__calc_sand_rho_sig_sand_num_1_5,
                  variance__calc_rho_sig_nom_sand_num_1_5,
                  sd__calc_sand_rho_sig_sand_num_1_5
                  ])

    tabd.add_row(['sand',
                  '6-10',
                  rmse__calc_sand_theta_nom_sand_num_6_10,
                  mb__calc_sand_theta_nom_sand_num_6_10,
                  pd__calc_sand_theta_nom_sand_num_6_10,
                  variance__calc_sand_theta_nom_sand_num_6_10,
                  sd__calc_sand_theta_nom_sand_num_6_10,

                  rmse__calc_sand_rho_nom_sand_num_6_10,
                  mb__calc_sand_rho_nom_sand_num_6_10,
                  pd__calc_sand_rho_nom_sand_num_6_10,
                  variance__calc_rho_theta_nom_sand_num_6_10,
                  sd__calc_sand_rho_nom_sand_num_6_10,

                  rmse__calc_sand_theta_sig_sand_num_6_10,
                  mb__calc_sand_theta_sig_sand_num_6_10,
                  pd__calc_sand_theta_sig_sand_num_6_10,
                  variance__calc_sand_theta_sig_sand_num_6_10,
                  sd__calc_sand_theta_sig_sand_num_6_10,

                  rmse__calc_sand_rho_sig_sand_num_6_10,
                  mb__calc_sand_rho_sig_sand_num_6_10,
                  pd__calc_sand_rho_sig_sand_num_6_10,
                  variance__calc_rho_sig_nom_sand_num_6_10,
                  sd__calc_sand_rho_sig_sand_num_6_10
                  ])

    tabd.add_row(['sand',
                  '11-15',
                  rmse__calc_sand_theta_nom_sand_num_11_15,
                  mb__calc_sand_theta_nom_sand_num_11_15,
                  pd__calc_sand_theta_nom_sand_num_11_15,
                  variance__calc_sand_theta_nom_sand_num_11_15,
                  sd__calc_sand_theta_nom_sand_num_11_15,

                  rmse__calc_sand_rho_nom_sand_num_11_15,
                  mb__calc_sand_rho_nom_sand_num_11_15,
                  pd__calc_sand_rho_nom_sand_num_11_15,
                  variance__calc_rho_theta_nom_sand_num_11_15,
                  sd__calc_sand_rho_nom_sand_num_11_15,

                  rmse__calc_sand_theta_sig_sand_num_11_15,
                  mb__calc_sand_theta_sig_sand_num_11_15,
                  pd__calc_sand_theta_sig_sand_num_11_15,
                  variance__calc_sand_theta_sig_sand_num_11_15,
                  sd__calc_sand_theta_sig_sand_num_11_15,

                  rmse__calc_sand_rho_sig_sand_num_11_15,
                  mb__calc_sand_rho_sig_sand_num_11_15,
                  pd__calc_sand_rho_sig_sand_num_11_15,
                  variance__calc_rho_sig_nom_sand_num_11_15,
                  sd__calc_sand_rho_sig_sand_num_11_15
                  ])

    tabd.add_row(['sand',
                  '16-20',
                  rmse__calc_sand_theta_nom_sand_num_16_20,
                  mb__calc_sand_theta_nom_sand_num_16_20,
                  pd__calc_sand_theta_nom_sand_num_16_20,
                  variance__calc_sand_theta_nom_sand_num_16_20,
                  sd__calc_sand_theta_nom_sand_num_16_20,

                  rmse__calc_sand_rho_nom_sand_num_16_20,
                  mb__calc_sand_rho_nom_sand_num_16_20,
                  pd__calc_sand_rho_nom_sand_num_16_20,
                  variance__calc_rho_theta_nom_sand_num_16_20,
                  sd__calc_sand_rho_nom_sand_num_16_20,

                  rmse__calc_sand_theta_sig_sand_num_16_20,
                  mb__calc_sand_theta_sig_sand_num_16_20,
                  pd__calc_sand_theta_sig_sand_num_16_20,
                  variance__calc_sand_theta_sig_sand_num_16_20,
                  sd__calc_sand_theta_sig_sand_num_16_20,

                  rmse__calc_sand_rho_sig_sand_num_16_20,
                  mb__calc_sand_rho_sig_sand_num_16_20,
                  pd__calc_sand_rho_sig_sand_num_16_20,
                  variance__calc_rho_sig_nom_sand_num_16_20,
                  sd__calc_sand_rho_sig_sand_num_16_20
                  ])

    tabd.add_row(['peat',
                  '21-25',

                  rmse__calc_peat_theta_nom_peat_num_21_25,
                  mb__calc_peat_theta_nom_peat_num_21_25,
                  pd__calc_peat_theta_nom_peat_num_21_25,
                  variance__calc_peat_theta_nom_peat_num_21_25,
                  sd__calc_peat_theta_nom_peat_num_21_25,

                  rmse__calc_peat_rho_nom_sand_num_21_25,
                  mb__calc_peat_rho_nom_sand_num_21_25,
                  pd__calc_peat_rho_nom_sand_num_21_25,
                  variance__calc_peat_rho_theta_nom_sand_num_21_25,
                  sd__calc_peat_rho_nom_sand_num_21_25,

                  rmse__calc_peat_theta_sig_sand_num_21_25,
                  mb__calc_peat_theta_sig_sand_num_21_25,
                  pd__calc_peat_theta_sig_sand_num_21_25,
                  variance__calc_peat_theta_sig_sand_num_21_25,
                  sd__calc_peat_theta_sig_sand_num_21_25,

                  rmse__calc_peat_rho_sig_peat_num_21_25,
                  mb__calc_peat_rho_sig_peat_num_21_25,
                  pd__calc_peat_rho_sig_peat_num_21_25,
                  variance__calc_peat_rho_sig_nom_peat_num_21_25,
                  sd__calc_peat_rho_sig_peat_num_21_25,
                  ])

    tabd.add_row(['peat',
                  '26-30',

                  rmse__calc_peat_theta_nom_peat_num_26_30,
                  mb__calc_peat_theta_nom_peat_num_26_30,
                  pd__calc_peat_theta_nom_peat_num_26_30,
                  variance__calc_peat_theta_nom_peat_num_26_30,
                  sd__calc_peat_theta_nom_peat_num_26_30,

                  rmse__calc_peat_rho_nom_peat_num_26_30,
                  mb__calc_peat_rho_nom_peat_num_26_30,
                  pd__calc_peat_rho_nom_peat_num_26_30,
                  variance__calc_peat_rho_theta_nom_peat_num_26_30,
                  sd__calc_peat_rho_nom_peat_num_26_30,

                  rmse__calc_peat_theta_sig_peat_num_26_30,
                  mb__calc_peat_theta_sig_peat_num_26_30,
                  pd__calc_peat_theta_sig_peat_num_26_30,
                  variance__calc_peat_theta_sig_peat_num_26_30,
                  sd__calc_peat_theta_sig_peat_num_26_30,

                  rmse__calc_peat_rho_sig_peat_num_26_30,
                  mb__calc_peat_rho_sig_peat_num_26_30,
                  pd__calc_peat_rho_sig_peat_num_26_30,
                  variance__calc_peat_rho_sig_nom_peat_num_26_30,
                  sd__calc_peat_rho_sig_peat_num_26_30
                  ])

    tabd.add_row(['peat',
                  '31-35',

                  rmse__calc_peat_theta_nom_peat_num_31_35,
                  mb__calc_peat_theta_nom_peat_num_31_35,
                  pd__calc_peat_theta_nom_peat_num_31_35,
                  variance__calc_peat_theta_nom_peat_num_31_35,
                  sd__calc_peat_theta_nom_peat_num_31_35,

                  rmse__calc_peat_rho_nom_peat_num_31_35,
                  mb__calc_peat_rho_nom_peat_num_31_35,
                  pd__calc_peat_rho_nom_peat_num_31_35,
                  variance__calc_peat_rho_theta_nom_peat_num_31_35,
                  sd__calc_peat_rho_nom_peat_num_31_35,

                  rmse__calc_peat_theta_sig_peat_num_31_35,
                  mb__calc_peat_theta_sig_peat_num_31_35,
                  pd__calc_peat_theta_sig_peat_num_31_35,
                  variance__calc_peat_theta_sig_peat_num_31_35,
                  sd__calc_peat_theta_sig_peat_num_31_35,

                  rmse__calc_peat_rho_sig_peat_num_31_35,
                  mb__calc_peat_rho_sig_peat_num_31_35,
                  pd__calc_peat_rho_sig_peat_num_31_35,
                  variance__calc__peat_rho_sig_nom_peat_num_31_35,
                  sd__calc_peat_rho_sig_peat_num_31_35
                  ])

    tabd.add_row(['SAND ALL',
                  '',
                  rmse__calc_sand_theta_nom_sand_all,
                  mb__calc_sand_theta_nom_sand_all,
                  pd__calc_sand_theta_nom_sand_all,
                  variance__calc_sand_theta_nom_sand_all,
                  sd__calc_sand_theta_nom_sand_trial_all,

                  rmse__calc_sand_rho_nom_sand_all,
                  mb__calc_sand_rho_nom_sand_all,
                  pd__calc_sand_rho_nom_sand_all,
                  variance__calc_sand_rho_nom_sand_all,
                  sd__calc_sand_rho_nom_sand_all,

                  rmse__calc_sand_theta_sig_sand_all,
                  mb__calc_sand_theta_sig_sand_all,
                  pd__calc_sand_theta_sig_sand_all,
                  variance__calc_sand_theta_sig_sand_all,
                  sd__calc_sand_theta_sig_sand_all,

                  rmse__calc_sand_rho_sig_sand_all,
                  mb__calc_sand_rho_sig_sand_all,
                  pd__calc_sand_rho_sig_sand_all,
                  variance__calc_sand_rho_sig_sand_all,
                  sd__calc_sand_rho_sig_sand_all
                  ])

    tabd.add_row(['PEAT ALL',
                  '',
                  rmse__calc_peat_theta_nom_peat_all,
                  mb__calc_peat_theta_nom_peat_all,
                  pd__calc_peat_theta_nom_peat_all,
                  variance__calc_peat_theta_nom_peat_all,
                  sd__calc_peat_theta_nom_peat_trial_all,

                  rmse__calc_peat_rho_nom_peat_all,
                  mb__calc_peat_rho_nom_peat_all,
                  pd__calc_peat_rho_nom_peat_all,
                  variance__calc_peat_rho_nom_peat_all,
                  sd__calc_peat_rho_nom_peat_all,

                  rmse__calc_peat_theta_sig_peat_all,
                  mb__calc_peat_theta_sig_peat_all,
                  pd__calc_peat_theta_sig_peat_all,
                  variance__calc_peat_theta_sig_peat_all,
                  sd__calc_peat_theta_sig_peat_all,

                  rmse__calc_peat_rho_sig_peat_all,
                  mb__calc_peat_rho_sig_peat_all,
                  pd__calc_peat_rho_sig_peat_all,
                  variance__calc_peat_rho_sig_peat_all,
                  sd__calc_peat_rho_sig_peat_all
                  ])

    print('Saving stat table...')
    html = tabd.get_html_string()
    fn = TABLE_PATH + TABLE_NAME_ALL + HTML_EXT
    with open(fn, 'w') as f:
        f.write(html)
        f.close()
    print('DONE saving out HTML table')

    # load in the HTML table text and transpose cols
    print('Converting to CSV')
    fn_csv = TABLE_PATH + TABLE_NAME_ALL + CSV_EXT
    df_list = pd.read_html(fn)
    for i, df in enumerate(df_list):
        df.to_csv(fn_csv)
    print('Done conversion to CSV')

    # load in the CSV file
    frame = np.genfromtxt(fn_csv, delimiter=',', dtype=str)
    frame1 = np.delete(frame, 0, 1)         # delete the first column
    frame2 = np.delete(frame1, 0, 0)        # delete the first row
    c0 = frame2[:, 0]
    c1 = frame2[:, 1]
    c2 = frame2[:, 2]
    c3 = frame2[:, 3]
    c4 = frame2[:, 4]
    c5 = frame2[:, 5]
    c6 = frame2[:, 6]
    c7 = frame2[:, 7]
    c8 = frame2[:, 8]
    c9 = frame2[:, 9]
    c10 = frame2[:, 10]
    c11 = frame2[:, 11]
    c12 = frame2[:, 12]
    c13 = frame2[:, 13]
    c14 = frame2[:, 14]
    c15 = frame2[:, 15]
    c16 = frame2[:, 16]
    c17 = frame2[:, 17]
    c18 = frame2[:, 18]
    c19 = frame2[:, 19]
    c20 = frame2[:, 20]
    c21 = frame2[:, 21]
    frame_out = np.column_stack((c0, c1,
                                 c2, c12,
                                 c3, c13,
                                 c4, c14,
                                 c5, c15,
                                 c6, c16,
                                 c7, c17,
                                 c8, c18,
                                 c9, c19,
                                 c10, c20,
                                 c11, c21))

    print('Saving sorted cols out as CSV file...')
    fn_csv_col_sorted = TABLE_PATH + TABLE_NAME_ALL + "-col-sorted" + CSV_EXT
    df = pd.DataFrame(frame_out)
    df.to_csv(fn_csv_col_sorted)

    fn_csv_col_sorted_html = TABLE_PATH + TABLE_NAME_ALL + "-col-sorted" + HTML_EXT
    print('Save to HTML file')
    df.to_html(fn_csv_col_sorted_html)
    print('DONE')


def main():
    data_analysis_for_trials()


if __name__ == '__main__':
    main()




