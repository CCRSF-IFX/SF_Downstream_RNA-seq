# Any Python script in the scripts folder will be able to import from this module.i
from snakemake.utils import validate
import pandas as pd
import os
import glob
configfile: "/config/config.yaml"
# Ensure the keys are available in the config
results = config['results']
raw_counts = config['raw_counts']
coldata = config['coldata']
