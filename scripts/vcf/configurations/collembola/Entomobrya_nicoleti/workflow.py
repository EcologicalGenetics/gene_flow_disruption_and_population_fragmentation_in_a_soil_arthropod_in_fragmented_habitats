#!/bin/env python3
import sys, os
sys.path.insert(0, os.path.realpath('../../../workflow_source/'))
from workflow_source import *

gwf = freebayes_population_set_workflow()
