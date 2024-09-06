#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Download Klebsiella data from Pathongen Watch
curl -L -o - https://pathogenwatch-public.ams3.cdn.digitaloceanspaces.com/Klebsiella%20pneumoniae__kleborate.csv.gz | gunzip > "${SCRIPT_DIR}/data/klebsiella_pneumoniae__kleborate.csv";
curl -L -o - https://pathogenwatch-public.ams3.cdn.digitaloceanspaces.com/Klebsiella%20pneumoniae__metadata.csv.gz | gunzip > "${SCRIPT_DIR}/data/klebsiella_pneumoniae__metadata.csv";