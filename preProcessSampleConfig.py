import os, sys
import pandas as pd

def baseNameFormat(df, cols, delim = "-"):
    # Input: dataframe, cols = colnames in order of basename
    # delim = string delimiter for basename column values
   
    # Checks:
    #
    # -- cols are members of df:
    nonExistColnames = [c for c in cols if c not in list(df)]
    if len(nonExistColnames) != 0:
        errorString = "Some columns in 'cols' don't exist:\n" + " ".join(nonExistColnames)
        raise ValueError(errorString)

    ###
    # make format string:
    baseNameFormat = delim.join(["{}" for c in cols])
    return(baseNameFormat)

def baseNameFromCols(df, cols, formatString, outputColName = "basename"):
    # Appends basename column to dataframe
    # Input:
    df[outputColName] = df[cols].apply(lambda x : formatString.format(*x), axis = 1)
    return(df)

def addBaseName(df, cols, delim = "-", outputColName = "basename"):
    # Takes sampleSheet & list of columns as input, 
    # returns sampleSheet with 'basename' column appended
    formatString = baseNameFormat(df, cols, delim)
    df = baseNameFromCols(df, cols, formatString, outputColName = outputColName)
    return(df)

#def addExt(df, ext, baseNameCol = 'basename'):
#    # input sampleSheet, return extension of file


def main(path, idcols, delim):
    df = pd.read_table(path, delimiter = "\t")
    df = addBaseName(df, idcols, delim)
    #return(df, pool_df)
    keep_cols = idcols.copy()
    keep_cols.append('basename')
    pool_df = df[keep_cols].drop_duplicates()
    
    return(df, pool_df)

# Drop call from argv
sys.argv.pop(0)

# Path to configFile is first argument
path = sys.argv.pop(0)
# TODO: check that path exists ?

# Remaining function calls are the id variables for the basename
idcols = sys.argv

sampleSheet, poolSampleSheet = main(path, idcols, '-')

out = {"sampleSheet": sampleSheet, "poolSampleSheet": poolSampleSheet}

for df, name in zip(out.values(), out.keys()):
    filename = name + ".tsv"
    df.to_csv(filename, sep = "\t", index = False)
