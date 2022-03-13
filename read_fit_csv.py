import numpy
import math
import scipy.io

radius = 2.105/(2*math.pi);
timestamps = []
omegas = []

with open('test.csv', newline='') as csvfile:
    lines = csvfile.readlines()
for line in lines:
    l = line.split(',')
    idx_timestamp = [i for i, x in enumerate(l) if x=='timestamp']
    idx_speed = [i for i, x in enumerate(l) if x=='speed']
    if len(idx_speed)>0:
        t = float(l[ idx_timestamp[0]+1 ].strip('"') ) + 631065600
        if t > 1000000000.0:
            timestamps.append(t)
            omegas.append(float(
                l[ idx_speed[0]+1 ].strip('"') )/radius)

wheelspeed = {
        'timestamp':timestamps,
        'omega':omegas}
scipy.io.savemat('wheelspeed_meas.mat',wheelspeed)
