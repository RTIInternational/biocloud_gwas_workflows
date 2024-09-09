import argparse
import boto3
import json
from datetime import datetime
import os
from zipfile import ZipFile

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument(
    '--working_dir',
    help='Directory for temporary files',
    type = str,
    required = True
)
parser.add_argument(
    '--repo_dir',
    help='Location of repository containing workflow',
    type = str,
    required = True
)
parser.add_argument(
    '--zip_include',
    help='Filters for repoository files to include in the repository zip file',
    type = str,
    nargs = '*',
    required = True
)
parser.add_argument(
    '--s3_access_key',
    help='Access key for S3 bucket',
    type = str,
    required = True
)
parser.add_argument(
    '--s3_secret_access_key',
    help='Secret access key for S3 bucket',
    type = str,
    required = True
)
parser.add_argument(
    '--target_dir',
    help='Directory to copy files to',
    type = str,
    required = True
)
parser.add_argument(
    '--downloaded_samples',
    help='File containing list of previously samples - these samples are ignored in the source bucket',
    type = str,
    required = True
)
parser.add_argument(
    '--download_limit',
    help='Max number of samples to download',
    type = int,
    default = 1000,
    required = False
)
parser.add_argument(
    '--samples_to_download',
    help='(Optional) List of files to download',
    type = str,
    required = False
)
args = parser.parse_args()

# Create working directory if doesn't exist
working_dir = args.working_dir if (args.working_dir[-1] == "/") else (args.working_dir + "/")
os.system("mkdir -p {}".format(working_dir))

# Zip repo
repo_name = os.path.basename(args.repo_dir)
with ZipFile("{}{}.zip".format(working_dir, repo_name), 'w', zipfile.ZIP_DEFLATED) as myzip:
    
    myzip.write('testtext.txt')

# Zip repo
try:
    cmd = [

    ]
    result = subprocess.run(
        cmd,
        text=True,
        check=True,
        capture_output=True
    )
except subprocess.CalledProcessError as e:
    result = e

session = boto3.Session(aws_access_key_id=args.s3_access_key, aws_secret_access_key=args.s3_secret_access_key)
omics = session.client('omics')

with open('/home/ngaddis/data/temp/wgs_qc.zip', 'rb') as f:
    wf_def = f.read()

with open('/home/ngaddis/git/biocloud_gwas_workflows/wgs_qc/wdl_v1.1/wgs_qc_wf_step_2_parameters.json') as f:
    wf_params = json.load(f)

request_id = 'wgs_qc_step_2_' + str(datetime.now().timestamp())
response = omics.create_workflow(
    name='wgs_qc_step_2',
    description='Step 2 of WGS QC',
    engine='WDL',
    definitionZip=wf_def,
    main='wgs_qc_wf_step_2.wdl',
    parameterTemplate=wf_params,
    storageCapacity=2000,
    requestId=request_id
)
