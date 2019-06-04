from ca import CellularSpace
from parsepop import get_ca
import os
import numpy as np
import matplotlib.pyplot as plt

NUM_OF_RUNS = 10
NUM_OF_ITER = 200
TEST_DIR = "test/"

def run_test(density_map, beta, gamma, runs = NUM_OF_RUNS, iters = NUM_OF_ITER, weight = 7, d_inf = 7):
    params = {'beta': beta,
              'gamma': gamma,
              'wieght': weight,
              'time_step': 1,
              'd_inf': d_inf
              }
    cell = (12,7)
    
    dir_name =  TEST_DIR + "b{}g{}/".format(beta,gamma)   
    os.makedirs(dir_name, exist_ok=True)
    
    cellular_space = CellularSpace(params, density_map = density_map)
    
    for i in range(runs):
        print("Run number ",i)
        cellular_space.reset()
        cellular_space.get_infection(cell = cell, n_infected = 1);
        lines = list()
        for iter in range(iters):
            sus,inf,rec = cellular_space.get_stats()
            lines.append("{},{},{}".format(sus,inf,rec))
            cellular_space.update()
        with open(dir_name + "{}.csv".format(i), 'w') as file:
            file.write("\n".join(lines));
            
def plot_average(dir_name):
    stats = np.zeros((NUM_OF_ITER,3))
    for i in range(NUM_OF_RUNS):
        ls = np.genfromtxt(TEST_DIR + dir_name + "/{}.csv".format(i),delimiter=',')
        stats += ls
    stats /= NUM_OF_RUNS
    plt.plot(stats[:,0],'g')
    plt.plot(stats[:,1],'r')
    plt.plot(stats[:,2],'b')
    plt.savefig(TEST_DIR + '/{}.png'.format(dir_name))
    plt.clf()

if __name__ == '__main__':
    #density_map = np.array(get_ca('500x500data.csv'))
    #density_map = np.flip(density_map, axis = 0)
    
    values = [0.2,0.4,0.6,0.8,1]
    
    for beta in values:
        for gamma in values:
            #print("Running test for beta={} and gamma={}:".format(beta,gamma))
            #run_test(density_map,beta,gamma)
            plot_average("b{}g{}".format(beta,gamma))