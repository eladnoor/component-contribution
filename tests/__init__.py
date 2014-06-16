# -*- coding: utf-8 -*-
import logging
import sys
from time import time

formatter = logging.Formatter('%(asctime)s %(filename)s %(lineno)s %(levelname)s  %(message)s')

stdout_handler = logging.StreamHandler(sys.stdout)
stderr_handler = logging.StreamHandler(sys.stderr)
file_handler = logging.FileHandler('/tmp/component_contribution_%f.log'%time())
stdout_handler.setFormatter(formatter)
stderr_handler.setFormatter(formatter)
file_handler.setFormatter(formatter)

logger = logging.getLogger('')
logger.addHandler(stdout_handler)
logger.addHandler(file_handler)
logger.setLevel(logging.INFO)

