# load_data.py
from utils_prm_est import *
import pandas as pd
from instructions import data_path
import numpy as np
import matplotlib.pyplot as plt
import os
# load mat
import scipy.io as sio



def import_data(temps, crates):
    # csv file
    # first raw is header
    # header: Temperature, Crate, NormCap, V
    path = data_path()
    data = pd.read_csv(path, header=0, sep="\t")

    # select data for the specified temperatures and C-rates
    data = data.loc[data['Temperature'].isin(temps)]
    data = data.loc[data['Crate'].isin(crates)]
    # convert data to dictionary
    dict = {}
    for temp in temps:
        dict[temp] = {}
        data_temp = data.loc[data['Temperature'] == temp]
        for crate in crates:
            dict[temp][crate] = {}
            data_crate = data_temp.loc[data_temp['Crate'] == crate]
            dict[temp][crate]['V'] = np.array(data_crate['V'])
            dict[temp][crate]['NormCap'] = np.array(data_crate['NormCap'])

    return dict

# set of functions to convert data from different formats to the format used 
# in the comparison with the simulation data
def convert_maccor(maccor_data, data_folder, plot=False):
    df = pd.read_excel(maccor_data, skiprows=0, header=0, index_col=False)
    # charge = False
    # if charge:
    #     df = df.loc[df['Md'] == "C"]
    # else:
    #     df = df.loc[df['Md'] == "D"]
    df = df.loc[df['Md'] == "D"]

    min_crate = 0.1
    selected_cycles = [7,13]
    temp = 298

    ref_current = np.min(df["Current [A]"])
    theor1C = ref_current*(1/min_crate)
    theor1C = 0.00876
    theorCap = theor1C

    # create csv file 
    # header: Temperature, Crate, NormCap, V
    # take name from maccor_data
    name = os.path.basename(maccor_data)
    name_csv = name.split(".")[0] + ".csv"
    with open(os.path.join(data_folder, name_csv), 'w') as f:
        # first row is header
        f.write("Temperature\tCrate\tNormCap\tV\n")
        for i in df["Cycle No. combined"].unique():
            df_cyc = df.loc[df["Cycle No. combined"] == i][1:]
            current_of_that_cycle = np.max(df_cyc["Current [A]"])
            C_rate = current_of_that_cycle/theor1C
            C_rate = round(C_rate, 2)
            if C_rate > 0.99:
                C_rate = round(C_rate, 0)
            if C_rate == 0.34:
                C_rate = 0.33
            if i in selected_cycles:
                for j in df_cyc.index:
                    f.write(str(temp) + "\t" +
                            str(C_rate) + "\t" +
                            str(df_cyc["Cap. [Ah]"][j]/theorCap) + "\t" +
                            str(df_cyc["Voltage [V]"][j])
                            + "\n")
    if plot:
        # plot from csv file just created
        # count how many different temperatures are in the file
        # header = ["Temperature", "Crate", "NormCap", "V"]
        data = pd.read_csv(os.path.join(data_folder, name_csv), header=0, sep="\t")
        temps = data["Temperature"].unique()
        fig, ax = plt.subplots(len(temps), 1, sharex=True, sharey=True)
        i = 0
        for temp in temps:
            data_loc = data.loc[data["Temperature"] == temp]
            crates = data_loc["Crate"].unique()
            for crate in crates:
                data_loc_crate = data_loc.loc[data_loc["Crate"] == crate]
                if len(temps) == 1:
                    ax.plot(data_loc_crate["NormCap"],
                            data_loc_crate["V"], label=crate)
                else:
                    ax[i].plot(data_loc_crate["NormCap"],
                               data_loc_crate["V"], label=crate)
            if len(temps) == 1:
                ax.legend()
                ax.set_title(str(temp) + " K")
            else:
                ax[i].legend()
                ax[i].set_title(str(temp) + " K")
            i += 1
        plt.show()

def convert_lanha_txt(lanha_folder, data_folder, plot=False):
    lanha_data = "file.txt"
    df = pd.read_csv(os.path.join(lanha_folder, lanha_data), sep="\t", header=2)
    # ignore last column
    df = df.loc[df.iloc[:,8] == "D"]
    min_crate = 0.2
    temp = 298
    # ref_current = -np.min(df.iloc[:,6])
    ref_current = -0.0000297391
    theor1C = ref_current*(1/min_crate)*0.941
    theorCap = 148.75*1e-6
    theor1C = -theorCap
    cycles = [3,5,7,9,11,13]

    with open(os.path.join(data_folder, "LTO.csv"), 'w') as f:
        f.write("Temperature\tCrate\tNormCap\tV\n")
        for i in df.iloc[:,0].unique():
            df_cyc = df.loc[df.iloc[:,0] == i]
            currents = df_cyc.iloc[:,6].unique()
            number_curr = np.size(currents)
            curr = currents[int(number_curr/2)]
            # current_of_that_cycle = (df_cyc.iloc[:,6])[int(np.size(df_cyc.iloc[:,6])/2)]
            C_rate = curr/theor1C
            C_rate = round(C_rate, 2)
            # if C_rate > 0.99:
            #     C_rate = round(C_rate, 0)
            # if C_rate == 0.34:
            #     C_rate = 0.33
            if i in cycles:
                for j in df_cyc.index:
                    f.write(str(temp) + "\t" +
                            str(C_rate) + "\t" +
                            str(df_cyc.iloc[:,4][j]/theorCap) + "\t" +
                            str(df_cyc.iloc[:,7][j])
                            + "\n")
    if plot:
        # plot from csv file just created
        # count how many different temperatures are in the file
        # header = ["Temperature", "Crate", "NormCap", "V"]
        data = pd.read_csv(os.path.join(data_folder, 'LTO.csv'), header=0, sep="\t")
        temps = data["Temperature"].unique()
        fig, ax = plt.subplots(len(temps), 1, sharex=True, sharey=True)
        i = 0
        for temp in temps:
            data_loc = data.loc[data["Temperature"] == temp]
            crates = data_loc["Crate"].unique()
            for crate in crates:
                data_loc_crate = data_loc.loc[data_loc["Crate"] == crate]
                if len(temps) == 1:
                    ax.plot(data_loc_crate["NormCap"],
                            data_loc_crate["V"], label=crate)
                else:
                    ax[i].plot(data_loc_crate["NormCap"],
                               data_loc_crate["V"], label=crate)
            if len(temps) == 1:
                ax.legend()
                ax.set_title(str(temp) + " K")
            else:
                ax[i].legend()
                ax[i].set_title(str(temp) + " K")
            i += 1
        plt.show()
        

    return 0

def convert_csv_figures(csv_folderm, data_folder, plot = True):
    df = pd.read_csv(csv_folderm, sep=",", header=0)
    ff = np.array([])
    for f in df['x']:
        # convert from , to . 
        f = f.replace(',','.')
        ff = np.append(ff, float(f))
    # take '1C' column and make a numpy array
    volt_1C = np.array([])
    for v in df['1C']:
        # convert from , to . 
        v = v.replace(',','.')
        volt_1C = np.append(volt_1C, float(v))
    index = np.where(volt_1C > 2)[0]
    ff_1C = ff[index]
    volt_1C = volt_1C[index]

    volt_2C = np.array([])
    for v in df['2C']:
        # convert from , to . 
        v = v.replace(',','.')
        volt_2C = np.append(volt_2C, float(v))
    index = np.where(volt_2C > 2)[0]
    volt_2C = volt_2C[index]
    ff_2C = ff[index]


    volt_5C = np.array([])
    for v in df['5C']:
        # convert from , to . 
        v = v.replace(',','.')
        volt_5C = np.append(volt_5C, float(v))
    index = np.where(volt_5C > 2)[0]
    volt_5C = volt_5C[index]
    ff_5C = ff[index]

    volt_10C = np.array([])
    for v in df['10C']:
        # convert from , to . 
        v = v.replace(',','.')
        volt_10C = np.append(volt_10C, float(v))
    index = np.where(volt_10C > 2)[0][:-21]
    print(index)
    volt_10C = volt_10C[index]
    ff_10C = ff[index]

    with open(os.path.join(data_folder, "LMFP_06.csv"), 'w') as f:
        f.write("Temperature\tCrate\tNormCap\tV\n")
        for i in range(len(ff_1C)):
            f.write(str(298) + "\t" +
                    str(1) + "\t" +
                    str(ff_1C[i]) + "\t" +
                    str(volt_1C[i])
                    + "\n")
        for i in range(len(ff_2C)):
            f.write(str(298) + "\t" +
                    str(2) + "\t" +
                    str(ff_2C[i]) + "\t" +
                    str(volt_2C[i])
                    + "\n")
        for i in range(len(ff_5C)):
            f.write(str(298) + "\t" +
                    str(5) + "\t" +
                    str(ff_5C[i]) + "\t" +
                    str(volt_5C[i])
                    + "\n")
        for i in range(len(ff_10C)):
            f.write(str(298) + "\t" +
                    str(10) + "\t" +
                    str(ff_10C[i]) + "\t" +
                    str(volt_10C[i])
                    + "\n")
            

    if plot:
        # plot from csv file just created
        # count how many different temperatures are in the file
        # header = ["Temperature", "Crate", "NormCap", "V"]
        data = pd.read_csv(os.path.join(data_folder, 'LMFP_06.csv'), header=0, sep="\t")
        temps = data["Temperature"].unique()
        fig, ax = plt.subplots(len(temps), 1, sharex=True, sharey=True)
        i = 0
        for temp in temps:
            data_loc = data.loc[data["Temperature"] == temp]
            crates = data_loc["Crate"].unique()
            for crate in crates:
                data_loc_crate = data_loc.loc[data_loc["Crate"] == crate]
                if len(temps) == 1:
                    ax.plot(data_loc_crate["NormCap"],
                            data_loc_crate["V"], label=crate)
                else:
                    ax[i].plot(data_loc_crate["NormCap"],
                               data_loc_crate["V"], label=crate)
            if len(temps) == 1:
                ax.legend()
                ax.set_title(str(temp) + " K")
            else:
                ax[i].legend()
                ax[i].set_title(str(temp) + " K")
            i += 1
        plt.show()


# csv_folder = '/Users/pierfrancescoombrini/Desktop/LMFP/Dynamics/figure_voltages/engineering/y=06.csv'
# data_folder = '/Users/pierfrancescoombrini/Desktop/LMFP/Dynamics/figure_voltages/engineering'
# convert_csv_figures(csv_folder, data_folder)
# lanha_folder = r"C:\Users\pierfrancescoo\Desktop"
# data_folder = r'C:\Users\pierfrancescoo\Documents\Phase-field\common-mpet\mpet-dev\LTO_mem\datacsv'
# convert_lanha_txt(lanha_folder, data_folder, plot=True)
# mat_data_fold = r"C:\Users\pierfrancescoo\Documents\Phase-field\Plot\plotter-TUD\Plotter"
# data_folder = r'C:\Users\pierfrancescoo\Documents\Phase-field\common-mpet\mpet-dev\LFP_mem\datacsv'
# convert_mat(mat_data_fold, data_folder)
# maccor_data = r'C:\Users\pierfrancescoo\Documents\Phase-field\Plot\plotter-TUD\Plotter\Mark_NMC_data\EPI54_new_fast_discharge.xlsx'
# data_folder = r'C:\Users\pierfrancescoo\Documents\Phase-field\common-mpet\mpet-dev\NMC_Mark\datacsv'
# convert_maccor(maccor_data, data_folder, plot=True)
