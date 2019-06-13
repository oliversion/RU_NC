#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from parsepop import get_ca
import ca
import argparse

def prob_float(string):
    value = float(string)
    if (value < 0) or (value > 1):
        msg = "%r is not between 0 and 1" % string
        raise argparse.ArgumentTypeError(msg)
    return value
        
def pos_int(string):
    value = int(string)
    if not (value > 0):
        msg = "%r is not a positive integer" % string
        raise argparse.ArgumentTypeError(msg)
    return value

parser = argparse.ArgumentParser(
            description='Runs simulation of infection spread.',
            epilog='''Some examples of usage:\n

simulate.py
       Runs simluation with default values.
        - beta = 0.6
        - gamma = 0.6
        - weight of same cell on infection = 7
        - minimum days of infection = 5
        - vaccination rate = 0%
        - starting cell: x = 6 y = 12 number of infected = 20
        - length of simulation: 200 days

simulate.py --map
       Shows the density map of the Netherlands. Can be used for choosing the
       starting cells for infection with -s.
       
simulate.py -b 0.8 -g 0.4 -i 3
       Uses beta=0.8 and gamma=0.4 with minimum of 3 days of infection for
       simulation.
       
simulate.py -d 500
       Simulation will run for 500 days with default values.
       
simulate.py -s 10 10 15 -s 6 12 40
       Starts simulation with two cells infected:
            - cell on x=10, y=10 with 15 infected people
            - cell on x=6, y=12 with 40 infected people

simulate.py -v 0.9
       Simulation starts with 90% vaccinated people (= as resistant).
''', formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-m','--map',action='store_true',help='Program will show only the density map, will not run a simulation of infection.')
parser.add_argument('-b','--beta',help='Beta parameter of simulation influencing the strength of infection spread.',default=0.6,type=prob_float)
parser.add_argument('-g','--gamma',help='Probability of infected individual to recover.',default=0.6,type=prob_float)
parser.add_argument('-d','--days',help='Number of days of simulation.',default=200,type=pos_int)
parser.add_argument('-v','--vrate',help='The rate of vaccination (number between 0 and 1) that determines the initial ratio of vaccinated population.',default=0,type=prob_float)
parser.add_argument('-i','--idays',help='Minimum number of days infected individual stays infected before they have chance to recover.',default=5,type=pos_int)
parser.add_argument('-s','--start',nargs=3,action='append',help='Adds starting cell of infection. Has 3 arguments: x and y coordinates of cell and number of infected people. Can be used multiple times.',type=pos_int, metavar=('X','Y','INF'))
parser.add_argument('-w','--weight',help='The weight that infection spread of current cell have compared to neighbouring cells.',default=7,type=pos_int)

args = parser.parse_args()

density_map = np.array(get_ca('100x100data.csv'))
density_map = np.flip(density_map, axis = 0)

params = {'beta': args.beta,
          'gamma': args.gamma,
          'wieght': args.weight,
          'time_step': 1,
          'd_inf': args.idays
          }

cellular_space = ca.CellularSpace(params, x=24, y=24, density_map = density_map, vaccination_rate = args.vrate)

fig, ax = plt.subplots()

if args.map:
    cellular_space.plot(mode = 'dm', ax = ax)
    plt.show()
    exit()
    
if args.start == None:
    args.start = [[6,12,20]]
    
for start_cell in args.start:
    cell = (start_cell[1],start_cell[0])
    cellular_space.get_infection(cell = cell, n_infected = start_cell[2])
    
for it in range(args.days):
    ax.cla()
    cellular_space.plot(mode = 'sm', it = it, ax = ax)
    cellular_space.print_stats(it = it)
    plt.pause(0.01)
    cellular_space.update()