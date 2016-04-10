#!/usr/bin/env python3

import tables as tb
import numpy as np

f=tb.open_file('./test.h5', 'r')
frames=f.root.test.simulation_frames
