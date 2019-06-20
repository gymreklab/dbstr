#!/bin/bash

# Installing
apt-get update
apt-get -y install python3 python3-pip wget git-core sqlite3 libsqlite3-dev
pip3 install argparse flask dash dash_core_components dash_html_components dash_table_experiments pandas plotly numpy pyfaidx


git clone https://github.com/gymreklab/dbstr
cd dbstr
git checkout release

# Upload data needed (using scp for now. Server from somewhere else?)
# /storage/resources/dbase/dbSTR/dbSTR.db
# /storage/resources/dbase/human/hg19/hg19.fa

# Run WebSTR
./WebSTR/WebSTR.py --host 0.0.0.0 --port 80
