import pandas as pd

# Function to import data, takes filepath and file name
def import_raw_data(path,file):
    # Load the CSV as a raw list of lines to find the row number
    with open(path+file, "r",encoding="ISO-8859-1") as f:
        for line_number, line in enumerate(f):
            if 'Cycle Nr.' in line:
                start_idx = line_number
                print(start_idx)

    # Read the CSV again starting from that row
    df = pd.read_csv(path+file, encoding="ISO-8859-1", skiprows=start_idx, skipfooter=4, engine='python')
    # Rename first three columns to remove special characters
    print(df.columns)
    # Errors won't be raised but should work for either degree symbol type
    df = df.rename(columns={"Cycle Nr.": "Cycle_nr", "Time [s]": "Time", "Temp. [°C]": "Temp", "Temp. [Â°C]": "Temp"})
    
    #df = df.rename(columns={"Cycle Nr.": "Cycle_nr", "Time [s]": "Time", "Temp. [Â°C]": "Temp"},errors="raise")
    return df

# Paths to raw_data and long_data folders, should not change
path = "/home/callum/rotation2/plate_reader/raw_data/"
outpath = "/home/callum/rotation2/plate_reader/long_data/"

# Files being processed
files = ["2025-05-08-Silwett-77-DC3000-B728a(Result sheet).csv",
        "2025-05-14-Silwett-77-1448a-B728a-DC3000.csv",
        "2025-06-12-ctrEDTA-DC3000.csv",
]
# Corresponding experiment names
exp_names = ["silwet01",
             "silwet02",
             "ctrEDTA02",
             ]

# Set index of interest to run through code
# E.g. for silwet01, index = 0
exp_index = 2

# For each experiment set up may have variation in column and row conditions
# So the following blocks set up lists to access conditions


df = import_raw_data(path, files[exp_index])
print(df.tail(n=5))
# print(list(df.columns[3:]))
melted_df = pd.melt(df, id_vars=['Cycle_nr','Time','Temp'],value_vars=list(df.columns[3:]))
melted_df["row"] = melted_df["variable"].apply(lambda x : x[0])
melted_df["column"] = melted_df["variable"].apply(lambda x: x[1:])
melted_df["experiment"] = exp_names[exp_index]

# Now need to set the mapping for the experiment

if exp_index in [2]:
    if exp_index == 2:
        species = 'DC3000'
    else:
        species = '1448a'
    def assign_species(row):
        if row['row'] in ['G', 'H']:
            return 'blank'
        else:
            global species
            return species
    melted_df['species'] = melted_df.apply(assign_species, axis=1)

if exp_index > 1:
    # Do ion concentration
    def assign_ion_conc(row):
        if int(row['column']) < 7:
            return 2
        else:
            return 0.4
    melted_df['ion_conc'] = melted_df.apply(assign_ion_conc, axis=1)
    # Do chelator concentration (same for all columns)
    initial_conc = 10
    ion_dict = {'A':initial_conc,
                'B':initial_conc/5,
                'C':initial_conc/(5**2),
                'D':initial_conc/(5**3),
                'E':initial_conc/(5**4),
                'F':0,
                'G':0,
                'H':initial_conc,}
    melted_df['chelator_conc'] = melted_df["row"].apply(lambda x: ion_dict[str(x)])
    # Do chelator type (citrate or EDTA)
    def assign_chelator_type(row):
        if int(row['column']) < 4:
            return "Citrate"
        elif int(row['column']) < 7:
            return "EDTA"
        elif int(row['column']) < 10:
            return "Citrate"
        elif int(row['column']) < 13:
            return "EDTA"
    melted_df['chelator_type'] = melted_df.apply(assign_chelator_type, axis=1)
if exp_index in [2]:
    melted_df['ion_type'] = "Mg"

# Save long form dataframe for ggplot analysis
pd.DataFrame.to_csv(melted_df,outpath+exp_names[exp_index]+".csv",index=False)

# Setting percentage of silwet based on row
silwet10_dict = {'A':10,
                'B':1,
                'C':0.1,
                'D':0.01,
                'E':0.001,
                'F':0.0001,
                'G':0.00001,
                'H':0}
row_dicts = [silwet10_dict]
# Setting bacterial species based on column
silwet01_col_dict = {"1":'DC3000',
                     "2":'DC3000',
                     "3":'DC3000',
                     "4":'DC3000',
                     "5":'B728a',
                     "6":'B728a',
                     "7":'B728a',
                     "8":'B728a',
                     "12":'blank',}
silwet02_col_dict = {"1":'1448a',
                     "2":'1448a',
                     "3":'1448a',
                     "4":'B728a',
                     "5":'B728a',
                     "6":'B728a',
                     "7":'DC3000',
                     "8":'DC3000',
                     "9":'DC3000',
                     "10":'blank',
                     "11":'blank',
                     "12":'blank',}
column_dicts = [silwet01_col_dict,silwet02_col_dict]
# Set mapping dictionaries
if exp_index < 3:
    # First three experiment are 10 fold dilution series
    condition_index = 0
row_dict = row_dicts[condition_index]
column_dict = column_dicts[exp_index]

# Write condition and species using these mappings
melted_df["silwet_prc"] = melted_df["row"].apply(lambda x: row_dict[x])
# Numeric dictionary key wasn't working great so switched to str
melted_df["species"] = melted_df["column"].apply(lambda x: column_dict[str(x)])

# # Save long form dataframe for ggplot analysis
# pd.DataFrame.to_csv(melted_df,outpath+exp_names[exp_index]+".csv",index=False)


# Function method deprecated, but if all future experiments are similar would be nicer
# def classify_species(x):
#     x = int(x)
#     if x <= 3:
#         return "1448a"
#     elif x <= 6:
#         return "B728a"
#     elif x <= 9:
#         return "DC3000"
#     else:
#         return "blank"
# melted_df["species"] = melted_df["column"].apply(classify_species)
# print(melted_df.head(n=5))

# Outlier removal is deprecated, moved to R to streamline workflow
# def remove_outlier(dataframe, well):
#     """Removes dataframe rows with outlier well data"""
#     dataframe = dataframe[dataframe["variable"] != well]
#     return dataframe

#outlier_removed = remove_outlier(melted_df,"E8")


#.DataFrame.to_csv(outlier_removed,outpath+exp_names[exp_index]+"_r.csv",index=False)
