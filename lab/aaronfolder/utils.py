import pandas as pd
import numpy as np

#### Processing Data ####

def process_data(df, fml):
    """
    Processes data to be used in Mastering Metrics 1.1

    Parameters:
    df: a pandas dataframe
    fml: female? either 1 or 0

    Returns:
    list: a list of processed data to be used in format_data()
    """
    # where computed values are stored
    data_list = []

    data_list.append(fml)

    # the two datasets
    data_set_health_insurance = df[(df["hi"] == 1) & (df["fml"] == fml)]
    data_set_no_health_insurance = df[(df["hi"] == 0) & (df["fml"] == fml)]

    # computing the values for the health row
    health_no_insurance = np.round(np.array(data_set_no_health_insurance["hlth"],dtype=np.float).mean(),decimals=2)
    health_insurance = np.round(np.array(data_set_health_insurance["hlth"],dtype=np.float).mean(),decimals=2)
    health_diff = np.round(health_insurance - health_no_insurance, decimals = 2)
    health_insurance_sd = np.round(np.std(np.array(data_set_health_insurance["hlth"])),decimals=2)
    health_no_insurance_sd = np.round(np.std(np.array(data_set_no_health_insurance["hlth"])),decimals=2)
    health_diff_se = np.round(data_set_no_health_insurance["hlth"].sem()-(data_set_health_insurance["hlth"]).sem(),decimals=2)

    data_list.append(health_insurance)
    data_list.append(health_insurance_sd)
    data_list.append(health_no_insurance)
    data_list.append(health_no_insurance_sd)
    data_list.append(health_diff)
    data_list.append(health_diff_se)

    for column in ["nwhite", "age", "yedu", "famsize", "empl", "inc"]:
        # computes the values for all of the other columns
        data_point_health_insurance = np.round(np.mean(data_set_health_insurance[column]),decimals=2)
        data_point_no_health_insurance = np.round(np.mean(data_set_no_health_insurance[column]),decimals=2)
        data_difference = np.round(data_point_health_insurance - data_point_no_health_insurance, decimals=2)
        data_difference_se = np.round(data_set_no_health_insurance[column].sem()-data_set_health_insurance[column].sem(), decimals=2)

        data_list.append(data_point_health_insurance)
        data_list.append(data_point_no_health_insurance)
        data_list.append(data_difference)
        data_list.append(data_difference_se)

    # computes sample size
    data_length_health_insurance = len(data_set_health_insurance)
    data_length_no_health_insurance = len(data_set_no_health_insurance)

    data_list.append(data_length_health_insurance)
    data_list.append(data_length_no_health_insurance)

    return data_list

#### Formatting Data ####

def format_data(listed):
    """
    Formats data like Table 1.1 of Mastering Metrics

    Parameters:
    listed: a python list-like of processed data (must be output of process_data)

    Returns:
    Text and a pandas dataframe: prints out 'Husbands' or 'Wives' according to the input,
        displays the resulting formatted dataframe, and prints out the
        Table explanation.

    """
    if listed[0] == 0:
        print("Husbands")
    else:
        print("Wives")
    data_frame = pd.DataFrame(data={' ': ["Health index", "Nonwhite","Age","Education","Family Size","Employed","Family Income","Sample Size"],
                   'Some HI (1)': ["{0} [{1}]".format(listed[1],listed[2]), listed[7],listed[11],listed[15],listed[19],listed[23],listed[27],listed[31]],
                  'No HI (0)': ["{0} [{1}]".format(listed[3],listed[4]), listed[8],listed[12],listed[16],listed[20],listed[24],listed[28],listed[32]],
                  'Difference (3)': ["{0} ({1})".format(listed[5],listed[6]), "{0} ({1})".format(listed[9],listed[10]),"{0} ({1})".format(listed[13],listed[14]),"{0} ({1})".format(listed[17],listed[18]),"{0} ({1})".format(listed[21],listed[22]),"{0} ({1})".format(listed[25],listed[26]),"{0} ({1})".format(listed[29],listed[30])," "]})
    display(data_frame)
    print("""Notes: This table reports average characteristics for insured and uninsured married couples in the
          2009 National Health Interview Survey (NHIS). Columns (1), (2), (4), and (5) show average characteristics of
          the group of individuals specified by the column heading. Columns (3) and (6) report the difference between
          the average characteristic for individuals with and without health insurance (HI).
          Standard deviations are in brackets; standard errorsare reported in parentheses.""")
