import argparse
import os
import sqlite3

import pandas as pd

SEPARATOR = "\t"

EXPECTED_FILE_EXTENSION = ".tsv"

DB_TABLE_NAME = 'mirna_targets'

DROP_COLUMN_MARKER = "_drop"
COL_NAME_DICT = {'mirtarbase': {'miRTarBase ID'                : DROP_COLUMN_MARKER+'1', 
                                'miRNA'                        : 'mirna',
                                'Species (miRNA)'              : DROP_COLUMN_MARKER+'3',
                                'Target Gene'                  : 'target_gene',
                                'Target Gene (Entrez Gene ID)' : DROP_COLUMN_MARKER+'5',
                                'Species (Target Gene)'        : DROP_COLUMN_MARKER+'6',
                                'Experiments'                  : 'experiments',
                                'Support Type'                 : 'support_type',
                                'References (PMID)'            : DROP_COLUMN_MARKER+'9'
                                },
                 'mirwalk_3UTR': {'miRNA'                      : 'mirna', 
                                'mRNA'                         : DROP_COLUMN_MARKER+'2',
                                'Genesymbol'                   : 'target_gene',
                                'binding_site'                 : DROP_COLUMN_MARKER+'4',
                                'binding_probability'          : 'confidence'
                                }
                } 

def main(data_path, db_path):
    """
    Main function of script.

    :param data_path: directory containing *.tsv source files (script will try to use every *.tsv file)
    :param db_path: full (including filename) destination path of result database
    """
    source_files = [filename for filename in os.listdir(data_path) if filename.endswith(EXPECTED_FILE_EXTENSION)]
    
    for filename in source_files:
        
        source_db_name = get_source_db_from_filename(filename)
        
        if source_db_name not in COL_NAME_DICT:
            raise InvalidNameError("%s is not recognized as miRNA target database name" % source_db_name)
        
        mirna_df = extract_data_frame_from_file(data_path, filename) 
        reformatted_df = reformat_data_frame(mirna_df, 
                                            COL_NAME_DICT[source_db_name], 
                                            DROP_COLUMN_MARKER)      
    
        create_db_from_data_frame(db_path, source_db_name, reformatted_df)


def get_source_db_from_filename(filename):
    """ 
    (WARNING!!!: FILE SHOULD BE NAMED AFTER SOURCE DB NAME, FOR EXAMPLE "mirtarbase.tsv") 
    Source DB names must be listed in COL_NAME_DICT
    :param filename: name of source file 
    :return: name of the source DB
    """
    return filename[:-len(EXPECTED_FILE_EXTENSION)]


def extract_data_frame_from_file(data_path, filename):
    """
    :param data_path: path of source file directory
    :param filename: source file name
    :return: pandas dataframe with file content
    """
    file_path = os.path.join(data_path, filename)
    return pd.read_csv(file_path, skiprows=0, sep=SEPARATOR)


def reformat_data_frame(data_frame, col_mapping_dict, drop_column_marker):
    """
    :param data_frame: pandas dataframe
    :param col_mapping_dict: map of column names for renaming/dropping
    :param drop_column_marker: prefix of column names to drop
    :return: reformatted dataframe
    """
    # rename columns
    data_frame.columns = [col_mapping_dict[col_name] for col_name in data_frame.columns]
    # drop redundant columns
    data_frame = data_frame[[c for c in data_frame.columns if not c.startswith(drop_column_marker)]]
    
    # if confidence is provided, sort asecending by confidence and keep only interactions with highest confidence
    if 'confidence' in col_mapping_dict:
        data_frame = data_frame.sort_values(by='confidence')
        cols = list(data_frame.columns)
        cols.remove('confidence')
        data_frame.drop_duplicates(cols, keep='last')
    else:    
        # drop duplicate interactions
        data_frame = data_frame.drop_duplicates()
    
    return data_frame


def create_db_from_data_frame(db_path, table_name, df):
    """
    Creates database from  whole content which was collected and reformatted as dataframe.

    :param db_path: full path of the created database
    :param table_name: name of the table
    :param df: content of the table
    """
    conn = sqlite3.connect(db_path)
    df.to_sql(table_name, conn, if_exists="replace")


def get_input():
    """
    Very simple CLI provider.
    :return: user input
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i")
    parser.add_argument("-o")
    return parser.parse_args()


if __name__ == '__main__':
    args = get_input()
    print("input dir {} output file {}".format(args.i, args.o))
    main(args.i, args.o)
