#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as pl
import dynamic_graph as dg
import dynamic_graph.sot_solver as ss
import time
import datetime as dt
import rospy
import csv
import tf

# SolverKine solver and relavant tasks
from dynamic_graph.sot.core import *


if __name__ == '__main__' :
    a = ss.SotSolver("equalitysolver")
    a.testsolver()
