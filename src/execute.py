import subprocess
import argparse
import pandas as pd
import toml
import socket
from azure.identity import DefaultAzureCredential
from azure.storage.filedatalake import DataLakeServiceClient
import argparse
import os
import pandas as pd
from dotenv import load_dotenv
import shutil

def build_cmd_args() -> argparse.Namespace:
    """Command-line interface builder. 
    
    To run execute.py, we will require two arguments, contract_ID and request_ID.
    Example: python execute.py --contract_ID VRMC_OPS --request_ID EF4DD65B-B09B-47DB-ABF4-2EB8D9B2A211
    """
    
    parser = argparse.ArgumentParser(description="Execute the whole normal optimisation work flow.")

    parser.add_argument('--contract_ID', action='store', type=str, required=True)
    parser.add_argument('--request_ID', action='store', type=str, required=True)

    args = parser.parse_args()

    return args


def read_config(file_path='./config.toml') -> dict:
    """Parses a config.toml file returning the result as dict. 
    
    We use this to point to the virtual environment python executable when run locally.
    """

    with open(file_path, 'r') as f:
        config = toml.load(f)

    return config

def load_env_var():
    """Load environment variables. 
    
    If on Batch, the .env file is located at /opt/miniconda/envs/rdna/.env. Otherwise, place
    the .env file in the root folder.
    """
        
    if os.path.exists('/opt/miniconda/envs/rdna/.env'): # Batch .env location
        load_dotenv('/opt/miniconda/envs/rdna/.env')
    else: # Local .env location. 
        load_dotenv()


def extract(contract_ID: str, request_ID: str):
    """Extract the data of the given run ID from ADLS, save it locally, and return the optimisation start datetime.
    
    The endpoint will point towards the DEV accounts if ran locally. If ran on the Batch, it will point towards PROD.
    When saving the data, we create a new local directory /runs/request_ID/data.
    The start datetime is retrieved from the data, and passed to the subprocesses.
    """

    disp_bar = f"{'-'*40}"
    print(disp_bar)
    print(f"{'Extracting data from ADLS...' : ^40}")

    # Authenticate to Azure
    default_credential = DefaultAzureCredential(exclude_shared_token_cache_credential=True) # https://stackoverflow.com/questions/67165101/azure-chainedtokencredential-fails-after-password-change

    # Data lake service client object
    account_url = os.getenv("ADLS_BATCH_DFS_URL")
    dl_service_client = DataLakeServiceClient(account_url, credential=default_credential)

    # Download files from ADLS to local disk
    container = "data"

    files = ["v_GetOptimiseInputItems.parquet", "v_GetOptimiseInspectionInputItems.parquet", "v_DomainValueDepot.parquet", 
             "v_GetOptimiseInputCrews.parquet", "v_GetDistTimeCacheForOptRun.parquet"] # defect list, inspection list, depot list, crew list, matrix data
    src_directory = "dmroads/optimisation/input/" + contract_ID + "/" + request_ID # Extracting from ADLS storage
    dst_directory = f"../runs/{request_ID}/data/" # Saving to local directory

    print(f"Storage account: {account_url} \nContainer: {container} \nDirectory: {src_directory}\n")
    print(f"Saving extracted files to local directory {dst_directory}\n")

    if os.path.exists(dst_directory):
        shutil.rmtree(dst_directory)
    
    os.makedirs(dst_directory)

    for file in files:
        src_path = src_directory + "/" + file
        dst_path = dst_directory + "/" + file

        print(f"Extracting {file}, and copying to local storage path.")
        file_client = dl_service_client.get_file_client(container, src_path) # File client object

        with open(file=dst_path, mode="wb") as local_file:
            data = file_client.download_file()
            data.readinto(local_file)


    crew_input = pd.read_parquet(os.path.join(dst_directory, files[3]))
    opt_start_date = crew_input["OptStartDate"][0].to_pydatetime() # Convert timestamp to datetime

    s_date_str = opt_start_date.date().strftime("%Y-%m-%d")
    s_time_str = "00:00" # Hard coded opt run start time

    print("Extraction complete.")
    print(disp_bar)

    return s_date_str, s_time_str


if __name__ == '__main__':
    # read config file
    config = read_config()
    
    python_exe = config['binary_paths']['python'] # Project specific Python 
    if python_exe == "INPUT_PATH_TO_PYTHON_EXECUTABLE":
        print("\nA virtual environment Python executable path has not been specified in config.toml.\n")
    julia_exe = "julia"
    
    # get cmd args
    args = build_cmd_args()
    contract_ID = args.contract_ID
    request_ID = args.request_ID

    # extract the GetOptimiseInputItems, GetOptimiseInputCrews, DomainValueDepot and the start datetime
    s_date_str, s_time_str = "2023-10-17", "00:00"

    print(f"{'OPTIMISATION DETAILS' : ^40}")
    print(f"Optimisation run on {s_date_str}. \nOptRunID: {request_ID} \nContractID: {contract_ID}")

    # Looping through each day in the horizon
    subprocess.run([python_exe, "preprocessing.py", "--request_ID", f"{request_ID}", "--s_date", f"{s_date_str}", '--s_time', f'{s_time_str}'])
    print('\nBooting Julia...')
    subprocess.run([julia_exe, "optimise.jl", "--request_ID", f"{request_ID}", "--s_date", f"{s_date_str}", '--s_time', f'{s_time_str}'])
    subprocess.run([python_exe, "plot_solution.py", "--request_ID", f"{request_ID}"]) #, "--s_date", f"{s_date_str}", '--s_time', f'{s_time_str}'])
    # subprocess.run([python_exe, "load.py", "--request_ID", f"{request_ID}", "--contract_ID", f"{contract_ID}"])    

    # Save subprocess outputs to separate txt file if on Batch
    if not socket.gethostname().startswith("DOW"):
        terminal_outputs = ["stdout.txt", "stderr.txt"]
        # Authenticate to Azure
        account_url = os.getenv("ADLS_BATCH_DFS_URL")
        default_credential = DefaultAzureCredential(exclude_shared_token_cache_credential=True)

        # Data lake service client object
        dl_service_client = DataLakeServiceClient(account_url, credential=default_credential)

        # Upload local files to ADLS
        container = "data"

        dst_directory = "dmroads/optimisation/output/" + contract_ID + "/" + request_ID + "/day1"
        file_system_client = dl_service_client.get_file_system_client(container)
        file_system_client.create_directory(dst_directory)
        
        for file in terminal_outputs:
            src_path = "../../" + file
            dst_path = dst_directory + "/" + file

            directory_client = file_system_client.get_directory_client(dst_directory)
            file_client = directory_client.create_file(file)
            with open(src_path, "rb")as f:
                raw_bytes = f.read()
                print(f"Uploading {file}.")
                file_client.upload_data(raw_bytes, overwrite=True)