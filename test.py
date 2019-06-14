from ca import CellularSpace
from parsepop import get_ca
import os
import numpy as np
import matplotlib.pyplot as plt

NUM_OF_RUNS = 10
NUM_OF_ITER = 200
TEST_DIR = "test/"
FIGURES_DIR = "figures/"

def run_test(density_map, beta, gamma, vacrate, runs = NUM_OF_RUNS, iters = NUM_OF_ITER, weight = 7, d_inf = 5):
    params = {'beta': beta,
              'gamma': gamma,
              'wieght': weight,
              'time_step': 1,
              'd_inf': d_inf
              }
    cell = (12,6)
    
    dir_name =  TEST_DIR + "v{}b{}g{}/".format(vacrate,beta,gamma)   
    os.makedirs(dir_name, exist_ok=True)
    
    cellular_space = CellularSpace(params, density_map = density_map, vaccination_rate = vacrate)
    
    for i in range(runs):
        print("Run number ",i)
        cellular_space.reset()
        cellular_space.get_infection(cell = cell, n_infected = 20);
        lines = list()
        for iter in range(iters):
            sus,inf,rec = cellular_space.get_stats()
            lines.append("{},{},{}".format(sus,inf,rec))
            cellular_space.update()
        with open(dir_name + "{}.csv".format(i), 'w') as file:
            file.write("\n".join(lines));
    
def load_data(dir_name):
    stats = np.zeros((NUM_OF_RUNS,NUM_OF_ITER,3))
    for i in range(NUM_OF_RUNS):
        stats[i] = np.genfromtxt(TEST_DIR + dir_name + "/{}.csv".format(i),delimiter=',')
    return stats
    
def average_data(stats):
    return stats.sum(axis=0) / NUM_OF_RUNS
    
def plot_average(stats, dir_name):
    max = np.max(stats,axis=0)
    min = np.min(stats,axis=0)
    stats_av = average_data(stats)
    
    new_dir_name =  FIGURES_DIR + "sirplots/"  
    os.makedirs(new_dir_name, exist_ok=True)
    
    plt.plot(stats_av[:,0],'g', label='suspected')
    plt.plot(min[:,0],'g--')
    plt.plot(max[:,0],'g--')
    plt.plot(stats_av[:,1],'r', label='infected')
    plt.plot(min[:,1],'r--')
    plt.plot(max[:,1],'r--')
    plt.plot(stats_av[:,2],'b', label='resistant')
    plt.plot(min[:,2],'b--')
    plt.plot(max[:,2],'b--')
    plt.title('#groups vs days')
    plt.legend()
    plt.savefig(new_dir_name + '/{}.png'.format(dir_name))
    plt.clf()
    
def plot_first_day_infected(stats, dir_name):
    stats_av = average_data(stats)
    suspected = stats_av[:,0]
    firstdayinf = [0]
    for i in range(1,NUM_OF_ITER):
        firstdayinf.append(suspected[i-1] - suspected[i])
    
    new_dir_name =  FIGURES_DIR + "firstdayinfected/"  
    os.makedirs(new_dir_name, exist_ok=True)
    
    plt.title('First time infected vs days')
    plt.plot(firstdayinf,'r')
    plt.savefig(new_dir_name + '/{}.png'.format(dir_name))
    plt.clf()
    
    
    

if __name__ == '__main__':
    #density_map = np.array(get_ca('500x500data.csv'))
    #density_map = np.flip(density_map, axis = 0)
        
    values = [0.2,0.4,0.6,0.8,1]
    valuesvac = [0,0.5,0.6,0.7,0.8,0.85,0.9,0.95]
    
    finalcases = dict()
    
    for vacrate in valuesvac:
        for beta in values:
            for gamma in values:
                print("Running test for {}% vaccination rate, beta={} and gamma={}:".format(vacrate,beta,gamma))
                #run_test(density_map,beta,gamma,vacrate)
                
                dir_name = "v{}b{}g{}".format(vacrate,beta,gamma)
                stats = load_data(dir_name)
                #plot_average(stats, dir_name)
                #plot_first_day_infected(stats,dir_name)
                
                stats_av = average_data(stats)
                finalcases[(vacrate,beta,gama)] = 