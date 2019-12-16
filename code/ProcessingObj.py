import os
from recomp_theta_rho import recomp_theta_rho_additional
from constants import *


class ProcessingObj:

    def __init__(self, path, start_string, additional_text, q, t_heat, t_total, end_trim_time=0,
                 r_nominal=None, theta_o=None, theta_m=None, add_time_after_max=None):
        self.path = path
        self.start_string = start_string
        self.additional_text = additional_text
        self.q = q
        self.t_heat = t_heat
        self.t_total = t_total
        self.end_trim_time = end_trim_time
        self.outputs = None
        self.r_nominal = r_nominal
        self.theta_o = theta_o
        self.theta_m = theta_m
        self.add_time_after_max = add_time_after_max

    def get_info(self):
        return self.path, self.start_string, self.additional_text, self.q, \
               self.t_heat, self.t_total, self.end_trim_time, self.outputs

    def run(self):
        trim_sp_begin = 0  # set this to zero for most of the processing
        self.outputs = recomp_theta_rho_additional(self.path, self.start_string, self.additional_text,
                                                   self.q, self.t_heat, self.t_total, self.end_trim_time, trim_sp_begin,
                                                   self.r_nominal, self.theta_o, self.theta_m, self.add_time_after_max)

    def get_outputs(self):
        return self.outputs


