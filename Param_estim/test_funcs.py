from utils import *
import numpy as np

temps = [298]
protocols = ['pitt']
c_rates_cc = [[0.5,1,2]]  # C-rates
c_rates_gitt = [[1]]  # C-rates
cont_volts = [[3.376,3.352,3.271,3.221]]  # V

data_dic = import_data(temps, protocols, c_rates_cc, c_rates_gitt, cont_volts)