import re
import string
from string import Template
import os
import hjson
import yaml

__author__ = ["Ramón Ángel Ruiz Fernández"]
__email__ = ["rruizfer@cern.ch"]

SS_config = hjson.load(open('config/config_SS.hjson'))


include: 'selection/Snakefile'
include: 'time_resolution/Snakefile'
include: 'epm/Snakefile'
include: 'scq_cal/Snakefile'




