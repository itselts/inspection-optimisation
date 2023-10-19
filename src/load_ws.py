import argparse
from azure.identity import DefaultAzureCredential
from azure.storage.filedatalake import DataLakeServiceClient
import os
from dotenv import load_dotenv
import socket

if socket.gethostname().startswith("DOW"):
    load_dotenv()
else:
    load_dotenv('/opt/miniconda/envs/rdna/.env')


def build_cmd_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--contract_ID', action='store', type=str, required=True)
    parser.add_argument('--request_ID', action='store', type=str, required=True)

    args = parser.parse_args()

    return args

args = build_cmd_args()
request_ID = args.request_ID
contract_ID = args.contract_ID

disp_bar = f"{'-'*40}"
print(disp_bar)
print(f"{'Uploading results to ADLS...' : ^40}")

#### DATA LOAD TO AZURE #####
# Authenticate to Azure
account_url = os.getenv("ADLS_BATCH_DFS_URL")
default_credential = DefaultAzureCredential(exclude_shared_token_cache_credential=True)

# Data lake service client object
dl_service_client = DataLakeServiceClient(account_url, credential=default_credential)

# Upload local files to ADLS
container = "data"

src_directory = f"../runs/{request_ID}/warm_start"
dst_directory = "dmroads/optimisation/output/" + contract_ID + "/" + request_ID + "/warm_start"

file_system_client = dl_service_client.get_file_system_client(container)
file_system_client.create_directory(dst_directory)

# Uploading unsolved map
file = "folium_map.html"
src_path = f"../runs/{request_ID}/outputs" + "/" + file
dst_path = dst_directory + "/" + file

directory_client = file_system_client.get_directory_client(dst_directory)
file_client = directory_client.create_file(file)
with open(src_path, "rb")as f:
    raw_bytes = f.read()
    print(f"Uploading {file}.")
    file_client.upload_data(raw_bytes, overwrite=True)

# Uploading results
files = ["t.csv", "optrunoutput.csv", "folium_ws.html"]

for file in files:
    src_path = src_directory + "/" + file
    dst_path = dst_directory + "/" + file

    directory_client = file_system_client.get_directory_client(dst_directory)
    file_client = directory_client.create_file(file)
    with open(src_path, "rb")as f:
        raw_bytes = f.read()
        print(f"Uploading {file}.")
        file_client.upload_data(raw_bytes, overwrite=True)