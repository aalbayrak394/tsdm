import glob
import ntpath
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
import librosa

def librosa_read_wav_files(wav_files):
    
    if not isinstance(wav_files, list):
        wav_files = [wav_files]
    
    return [librosa.load(f)[0] for f in wav_files]


def data_generation():
    # List the wav files
    ROOT_DIR_TEST = glob.glob('./data/cats_dogs/test')[0]
    ROOT_DIR_TRAIN = glob.glob('./data/cats_dogs/train')[0]


    X_path = glob.glob(ROOT_DIR_TEST + "/test/*") # test = dogs in this case ! (wrong name of directory was given when it was created)
    X_path = X_path + glob.glob(ROOT_DIR_TEST + "/cats/*")
    X_path = X_path + glob.glob(ROOT_DIR_TRAIN + "/dog/*")
    X_path = X_path + glob.glob(ROOT_DIR_TRAIN + "/cat/*")

    y = np.empty((0, 1, ))
    for f in X_path:
        if 'cat' in ntpath.basename(f):
            resp = np.array([0])
            resp = resp.reshape(1, 1, )
            y = np.vstack((y, resp))
        elif 'dog' in ntpath.basename(f):
            resp = np.array([1])
            resp = resp.reshape(1, 1, )
            y = np.vstack((y, resp))
    
    
    print(f"in the data, there is {int(len(y) - sum(y))} cats (0) and {int(sum(y))} dogs (1)")

    wav_rate = librosa.load(X_path[0])[1]
    X = librosa_read_wav_files(X_path)
    
    return X, y, wav_rate


def get_data():
    
    df1_file_name = "covid.csv"
    df2_file_name = "Nowcasting_Zahlen.xlsx"
    
    df1 = pd.read_csv(df1_file_name)
    df2 = pd.read_excel(df2_file_name, sheet_name=1, engine='openpyxl')
    
    loc_df1 = df1.rename(columns={'dateRep': 'Date'})
    loc_df1['Date'] = pd.to_datetime(loc_df1.Date.values)
    loc_df1 = loc_df1.sort_values(by=['Date'])
    
    loc_df2 = df2.rename(columns={'Datum des Erkrankungsbeginns': 'Date'})
    loc_df2 = loc_df2.sort_values(by=['Date'])
    
    return loc_df1.set_index('Date'), loc_df2.set_index('Date')


def merge_two_datasets(df1, df2):

    merged_df = pd.merge_asof(df1, df2, on='Date').dropna()
    merged_df = merged_df.set_index('Date')
	
    return merged_df.sort_values(by=["year", "month", "day"])