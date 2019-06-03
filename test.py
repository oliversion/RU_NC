import ca
from parsepop import get_ca
import os

NUM_OF_RUNS = 10
NUM_OF_ITER = 10
TEST_DIR = "test/"

def run_test(density_map, beta, gamma, runs = NUM_OF_RUNS, iters = NUM_OF_ITER, weight = 7, d_inf = 5):
    params = {'beta': beta,
              'gamma': gamma,
              'wieght': weight,
              'time_step': 1,
              'd_inf': d_inf
              }
    cell = (12,7)
    
    dir_name =  TEST_DIR + "b{}g{}/".format(beta,gamma)   
    os.makedirs(dir_name)
    
    cellular_space = CellularSpace(params, density_map = density_map)
    
    for i in range(runs):
        cellular_space.reset()
        cellular_space.get_infection(cell = cell, n_infected = 1)
        lines = list()
        for iter in range(iters):
            s,i,r = cellular_space.get_stats()
            lines.append("{},{},{}".format(s,i,r))
            cellular_space.update()
        with open(dir_name + "{}.csv".format(i), 'w') as file:
            file.write("\n".join(lines));
            
        

if __name__ == '__main__':
    density_map = np.array(get_ca('500x500data.csv'))
    density_map = np.flip(density_map, axis = 0)
    
    values = [0.2,0.4,0.6,0.8,1]
    
    for beta in values:
        for gamma in values:
            run_test(density_map)