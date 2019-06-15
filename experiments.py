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
    
def plot_average(stats):
    max = np.max(stats,axis=0)
    min = np.min(stats,axis=0)
    stats_av = average_data(stats)
    vac_pop = stats_av[0,-1]
    
    plt.plot(stats_av[:,0] + vac_pop,'g', label='suspected')
    plt.plot(min[:,0] + vac_pop,'g--')
    plt.plot(max[:,0] + vac_pop,'g--')
    plt.plot(stats_av[:,1],'r', label='infected')
    plt.plot(min[:,1],'r--')
    plt.plot(max[:,1],'r--')
    plt.plot(stats_av[:,2] - vac_pop,'b', label='resistant')
    plt.plot(min[:,2] - vac_pop,'b--')
    plt.plot(max[:,2] - vac_pop,'b--')
    if vac_pop > 0:
        plt.plot(np.ones(NUM_OF_ITER)*vac_pop,'k', label='vaccinated')
    #plt.title('#groups vs days')
    plt.xlabel("days")
    plt.ylabel("population")
    plt.legend()
    
def plot_first_day_infected(stats, type='r', label='firstdayinf'):
    stats_av = average_data(stats)
    suspected = stats_av[:,0]
    firstdayinf = [0]
    for i in range(1,NUM_OF_ITER):
        firstdayinf.append(suspected[i-1] - suspected[i])
        
    #plt.title('First time infected vs days')
    plt.plot(firstdayinf,type,label = label)
    plt.xlabel("days")
    plt.ylabel("size of first time infected")
    
def plot_first_day_infected_rate(stats, type='r', label='firstdayinf'):
    stats_av = average_data(stats)
    suspected = stats_av[:,0]
    unvac_pop = stats_av[0,0]
    firstdayinf = [0]
    for i in range(1,NUM_OF_ITER):
        firstdayinf.append((suspected[i-1] - suspected[i])/unvac_pop)
        
    #plt.title('First time infected ratio against unvaccinated population vs days')
    plt.plot(firstdayinf,type,label = label)
    plt.xlabel("days")
    plt.ylabel(r'$R_I$',rotation=0)
    
def plot_num_cases_ratio(stats, type='r', label='num of cases'):
    stats_av = average_data(stats)
    suspected = stats_av[:,0]
    unvac_pop = stats_av[0,0]
    

    cases_ratio = (-suspected + unvac_pop)/unvac_pop
        
    #plt.title('Ratio of all time cases against unvaccinated population vs days')
    plt.plot(cases_ratio,type,label = label)
    plt.xlabel("days")
    plt.ylabel(r'$R_A$',rotation=0)
    
def save_fig(dir_name):
    plt.savefig(dir_name)
    plt.clf()
    
def mk_dir_and_save(dir_name,file_name):
    os.makedirs(dir_name, exist_ok=True)
    save_fig(dir_name + '/' + file_name)
    
def load_stats(beta,gamma,vacrate):
    file_name = "v{}b{}g{}".format(vacrate,beta,gamma)
    return load_data(file_name)
    
def plot_for_all(plot_function, dir_name):
    values = [0.2,0.4,0.6,0.8,1]
    valuesvac = [0,0.5,0.6,0.7,0.8,0.85,0.9,0.95]
    
    for vacrate in valuesvac:
        for beta in values:
            for gamma in values:
                stats = load_stats(beta,gamma,vacrate)
                
                plot_function(stats)
                mk_dir_and_save(FIGURES_DIR + dir_name, "v{}b{}g{}.png".format(vacrate,beta,gamma))
                
def plot_for_vac(plot_function, dir_name):
    values = [0.2,0.4,0.6,0.8,1]
    valsvactype = {0:'r', 0.5:'b', 0.8:'g', 0.85:'c', 0.9:'m', 0.95:'k'}
    
    for beta in values:
        for gamma in values:
            for vacrate in valsvactype:
                stats = load_stats(beta,gamma,vacrate)
                plot_function(stats,type=valsvactype[vacrate],label='vr = {}%'.format(vacrate*100))
            plt.legend()
            mk_dir_and_save(FIGURES_DIR + dir_name, "b{}g{}.png".format(beta,gamma))
            
def run_experiments():
    density_map = np.array(get_ca('500x500data.csv'))
    density_map = np.flip(density_map, axis = 0)
        
    values = [0.2,0.4,0.6,0.8,1]
    valuesvac = [0,0.5,0.6,0.7,0.8,0.85,0.9,0.95]
    
    for vacrate in valuesvac:
        for beta in values:
            for gamma in values:
                print("Running test for {}% vaccination rate, beta={} and gamma={}".format(vacrate,beta,gamma))
                run_test(density_map,beta,gamma,vacrate)
    
def get_reach(stats):
    stats_av = average_data(stats)
    starting_pop = stats_av[0,0]
    ending_pop = stats_av[-1,0]
    return (starting_pop - ending_pop)/starting_pop
    
def plot_reach(vacrate=0):
    values = [0.2,0.4,0.6,0.8,1]
    
    for gamma in values:
        x = list()
        y = list()
        for beta in values:
            stats = load_stats(beta,gamma,vacrate)
            x.append(beta)
            y.append(get_reach(stats))
        plt.plot(x,y,label=r'$\gamma$ = {}'.format(gamma))
    plt.legend()
    plt.xlabel(r'$\beta$')
    plt.ylabel('ratio of all time infected at the end')
    mk_dir_and_save(FIGURES_DIR, "reach.png")
    
if __name__ == '__main__':
    print('Running experiments')
    run_experiments()
    
    print('Creating plots')
    plot_for_all(plot_average,'sirplots')
    plot_for_all(plot_first_day_infected,'firstdayinfected')
    
    plot_for_vac(plot_first_day_infected_rate, "firstdaycomparison")
    plot_for_vac(plot_num_cases_ratio, "casesratio")
    
    plot_reach()