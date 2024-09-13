import argparse
import boto3
import json
from datetime import datetime
import os
from pathlib import Path
import zipfile

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument(
    '--working_dir',
    help = 'Directory for temporary files',
    type = str,
    required = True
)
parser.add_argument(
    '--repo_dir',
    help = 'Location of repository containing workflow',
    type = str,
    required = True
)
parser.add_argument(
    '--zip_include',
    help = 'Filters for repoository files to include in the repository zip file',
    type = str,
    required = True
)
parser.add_argument(
    '--s3_access_key',
    help = 'Access key for S3 bucket',
    type = str,
    required = True
)
parser.add_argument(
    '--s3_secret_access_key',
    help = 'Secret access key for S3 bucket',
    type = str,
    required = True
)
parser.add_argument(
    '--wf_parameters',
    help = 'JSON containing workflow parameters',
    type = str,
    required = True
)
parser.add_argument(
    '--wf_main_wdl',
    help = 'Main wdl file for workflow',
    type = str,
    required = True
)
parser.add_argument(
    '--wf_name',
    help = 'Name of wf to be created',
    type = str,
    required = True
)
parser.add_argument(
    '--wf_description',
    help = 'Description of wf to be created',
    type = str,
    required = True
)
parser.add_argument(
    '--wf_engine',
    help = 'Engine to use for workflow',
    type = str,
    default = 'WDL',
    choices = ['WDL', 'NEXTFLOW', 'CWL'],
    required = False
)
parser.add_argument(
    '--wf_storage_capacity',
    help = 'Storage capacity for workflow',
    type = int,
    default = 2000,
    choices = range(100001),
    required = False
)

args = parser.parse_args()

# Create working directory if doesn't exist
working_dir = args.working_dir if (args.working_dir[-1] == "/") else (args.working_dir + "/")
os.system("mkdir -p {}".format(working_dir))

# Zip repo
repo_name = os.path.basename(args.repo_dir)
zip_file = "{}{}.zip".format(working_dir, repo_name)
wdl_files = [str(p) for p in list(Path(args.repo_dir).rglob("*.wdl"))]
with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as myzip:
    for file in wdl_files:
        if args.zip_include in file:
            myzip.write(file, os.path.basename(file))

# Create Healthomics session
session = boto3.Session(aws_access_key_id=args.s3_access_key, aws_secret_access_key=args.s3_secret_access_key)
omics = session.client('omics')

with open(zip_file, 'rb') as f:
    wf_def = f.read()

with open(args.wf_parameters) as f:
    wf_params = json.load(f)

request_id = args.wf_name + str(datetime.now().timestamp())
response = omics.create_workflow(
    name=args.wf_name,
    description=args.wf_description,
    engine=args.wf_engine,
    definitionZip=wf_def,
    main=args.wf_main_wdl,
    parameterTemplate=wf_params,
    storageCapacity=args.wf_storage_capacity,
    requestId=request_id
)
