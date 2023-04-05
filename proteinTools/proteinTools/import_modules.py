import os
import xml.etree.ElementTree as ET
import sys
import csv
import mygene 
import requests
from functools import cached_property, lru_cache
import urllib
import pandas as pd
import numpy as np
import traceback
from io import StringIO
from chembl_webresource_client.new_client import new_client
import re