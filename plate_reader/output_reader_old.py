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
    df = pd.read_csv(path+file, encoding="ISO-8859-1", skiprows=start_idx, skipfooter=5, engine='python')
    
    # Rename first three columns to remove special characters
    print(df.columns)
    
    # Transpose
    df_t = df.T
    # Move the index into a column
    df_t = df_t.reset_index()
    
    # Set the top row as the column names
    # 1. Set the first row as the header
    df_t.columns = df_t.iloc[0]
    # 2. Drop the first row since it's now the header
    df_t = df_t[1:]
    # 3. (Optional) Reset the index if needed
    df_t = df_t.reset_index(drop=True)
        
    # Errors won't be raised but should work for either degree symbol type
    df_t = df_t.rename(columns={"Cycle Nr.": "Cycle_nr", "Time [s]": "Time", "Temp. [°C]": "Temp", "Temp. [Â°C]": "Temp"})
    print(df_t.head(5))


    return df_t

# Paths to raw_data and long_data folders, should not change
path = "/home/callum/rotation2/plate_reader/raw_data/"
outpath = "/home/callum/rotation2/plate_reader/long_data/"

# Files being processed
files = ["2025-05-08-Silwett-77-DC3000-B728a(Result sheet).csv",
        "2025-05-14-Silwett-77-1448a-B728a-DC3000.csv",
        "2025-05-22-slwtion-b728a-dc3000.csv",
        "2025-05-23-divion-b728a-dc3000.csv",
        "2025-05-30-divion-b728a-dc3000.csv",
        "2025-06-01-divion-b728a-dc3000.csv",
        "2025-06-04-ctrion-b728a-dc3000.csv",
        "2025-06-06-divion-b728a-dc3000.csv",
        "2025-06-11-ctrEDTA-b728a.csv",
        "2025-06-13-ctrEDTA-1448a.csv",
]
# Corresponding experiment names
exp_names = ["silwet01",
             "silwet02",
             "slwtion01",
             "divion01",
             "divion02",
             "divion03",
             "ctrion01",
             "divion04",
             "ctrEDTA01",
             "ctrEDTA03",
             ]

# Set index of interest to run through code
# E.g. for silwet01, index = 0
exp_index = 8

# For each experiment set up may have variation in column and row conditions
# So the following blocks set up lists to access conditions


df = import_raw_data(path, files[exp_index])
# print(df.tail(n=5))
# print(list(df.columns[3:]))
melted_df = pd.melt(df, id_vars=['Cycle_nr','Time','Temp'],value_vars=list(df.columns[3:]))
melted_df = melted_df.rename(columns = {0:'variable'})
melted_df["row"] = melted_df["variable"].apply(lambda x : x[0])
melted_df["column"] = melted_df["variable"].apply(lambda x: x[1:])
melted_df["experiment"] = exp_names[exp_index]

# Now need to set the mapping for the experiment
if exp_index < 8:
    def assign_species(row):
        if row['row'] in ['G', 'H']:
            return 'blank'
        elif int(row['column']) <= 6:
            return 'B728a'
        else:
            return 'DC3000'
    melted_df['species'] = melted_df.apply(assign_species, axis=1)
elif exp_index in [8,9]:
    if exp_index == 8:
        species = 'B728a'
    else:
        species = '1448a'
    def assign_species(row):
        if row['row'] in ['G', 'H']:
            return 'blank'
        else:
            global species
            return species
    melted_df['species'] = melted_df.apply(assign_species, axis=1)

initial_conc = 10

if exp_index in [4,5,7]:
    ion_dict = {'A':initial_conc,
                'B':initial_conc/5,
                'C':initial_conc/(5**2),
                'D':initial_conc/(5**3),
                'E':initial_conc/(5**4),
                'F':0,
                'G':0,
                'H':initial_conc,}
    melted_df['ion_conc'] = melted_df["row"].apply(lambda x: ion_dict[str(x)])
elif exp_index in [6]:
    ion_dict = {'A':initial_conc,
                'B':initial_conc/5,
                'C':initial_conc/(5**2),
                'D':initial_conc/(5**3),
                'E':initial_conc/(5**4),
                'F':0,
                'G':initial_conc,
                'H':0,}
    melted_df['citrate_conc'] = melted_df["row"].apply(lambda x: ion_dict[str(x)])
    def assign_ion_conc(row):
        if int(row['column']) < 4:
            return 2
        elif int(row['column']) < 7:
            return 0.4
        elif int(row['column']) < 10:
            return 2
        elif int(row['column']) < 13:
            return 0.4
    melted_df['ion_conc'] = melted_df.apply(assign_ion_conc, axis=1)
elif exp_index > 7:
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


ion_type = {'1':"Ca",
            '2':"Ca",
            '3':"Ca",
            '4':"Mg",
            '5':"Mg",
            '6':"Mg",
            '7':"Ca",
            '8':"Ca",
            '9':"Ca",
            '10':"Mg",
            '11':"Mg",
            '12':"Mg",}
if exp_index in [4,5,7]:
    melted_df['ion_type'] = melted_df["column"].apply(lambda x: ion_type[str(x)])
elif exp_index in [6,8,9]:
    melted_df['ion_type'] = "Mg"
# # Setting percentage of silwet based on column
# silwet10_dict = {'1':0,
#                 '2':0,
#                 '3':0,
#                 '4':10,
#                 '5':10,
#                 '6':10,
#                 '7':1,
#                 '8':1,
#                 '9':1,
#                 '10':0.1,
#                 '11':0.1,
#                 '12':0.1,}
# column_dicts = [silwet10_dict]
# # Setting bacterial species based on row
# species_dict = {"A":'B728a',
#                 "B":'DC3000',
#                 "C":'blank',
#                 "D":'B728a',
#                 "E":'DC3000',
#                 "F":'blank',
# }
# # Setting ion concentration based on row
# ion_dict = {"A":2,
#             "B":2,
#             "C":2,
#             "D":0.01,
#             "E":0.01,
#             "F":0.01,
# }


# # Write conditions and species using these mappings
# # Numeric dictionary key wasn't working great so switched to str
# melted_df["silwet_prc"] = melted_df["column"].apply(lambda x: silwet10_dict[str(x)])
# melted_df["species"] = melted_df["row"].apply(lambda x: species_dict[str(x)])
# melted_df["ion_mM"] = melted_df["row"].apply(lambda x: ion_dict[str(x)])

# Save long form dataframe for ggplot analysis
pd.DataFrame.to_csv(melted_df,outpath+exp_names[exp_index]+".csv",index=False)

