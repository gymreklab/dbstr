#!/bin/bash

# Mounting storage
# https://objectivefs.com/howto/how-to-mount-instance-store-on-aws-ec2-instance-for-disk-cache
sudo mkfs.ext4 -E nodiscard /dev/nvme1n1
sudo mkdir /var/cache/objectivefs
sudo mount /dev/nvme1n1 /var/cache/objectivefs

# Installing
sudo apt-get update
sudo apt-get -y install python3 python3-pip wget git-core sqlite3 libsqlite3-dev
sudo pip3 install argparse flask dash dash_core_components dash_html_components dash_table_experiments pandas plotly numpy pyfaidx

cd
git clone https://github.com/gymreklab/dbstr
cd dbstr
git checkout release

# Upload data needed (using scp for now. Server from somewhere else?)
# /storage/resources/dbase/dbSTR/dbSTR.db
# /storage/resources/dbase/human/hg19/hg19.fa

# Run WebSTR
sudo ./WebSTR/WebSTR.py --host 0.0.0.0 --port 80
