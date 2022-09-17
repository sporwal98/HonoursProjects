from program import *
def singlesite(output = 'output.txt', ss = 1e0, niterations = 1000, skip = 0):
    '''
    Main simulation code
    '''
    np.random.seed(42) #random seed for consistent results
    idt = 0 #euclidean time i in {1,...,ntau}
    tmc = 0 #Monte Carlo timesteps

    niter = niterations #number of iterations
    stepsize = ss #increments in position(permitted positions)
    x = 0.0 #initial position
    tmc = 0 #monte carlo time-step
        
    xs = [x] #list of x positions
    ts = [tmc] #list of time-steps
    
    '''
    #open output file
    outfile = open('output.txt','w')
    outfile.write('HEADERS')
    '''
    
    naccept = 0 #counter for acceptances
    sum_x = 0.0 #sum of x, used for <x>
    sum_xx = 0.0 #sum of x^2, used for <x^2>

    act = Action(m = 1, om = 1) #mass and omega just for generalisation
    met = Metropolis(ss = stepsize, act = act) 

    energies = [act.harmonic_act(x)]
    
    for iter in range(niter):
    #For each monte carlo timestep
        x, acc = met.step(x)
        tmc += 1
        

        if(acc):
            naccept += 1

        xs.append(x)
        ts.append(tmc)
        energies.append(act.harmonic_act(x))
        
        sum_x += x
        sum_xx += x**2
        #outfile.write('DATA')

    #outfile.close()

    '''   
    plt.clf()
    plt.plot(xs[skip:],ts[skip:])
    plt.xlabel('Position x')
    plt.ylabel('Monte Carlo Timestep')
    plt.title('Acceptance = '+str(naccept/niter))
    plt.show()

    plt.clf()
    plt.plot(ts[skip:], energies[skip:])
    plt.xlabel('Monte Carlo Timestep')
    plt.ylabel('Action')
    plt.title('Energies')
    plt.show()
    '''

    mu = sum_x/niter
    rms = sum_xx/niter
    var = rms - (mu**2)

    sigma = np.sqrt(var)
    gaussx = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
    normal = np.exp(-gaussx**2/2*var)/(sigma*np.sqrt(2*np.pi))
    
    
    plt.clf()
    plt.plot(gaussx, normal,color = 'r',label = 'gaussian')
    plt.hist(xs, bins = int(np.sqrt(niter)),label = 'histogram', density = True)
    plt.title('Position Histogram and Gaussian '+ str(niter))
    plt.legend()
    plt.show()
    
    print('<x> = ' + str(mu))
    print('<x^2> = ' + str(rms))
    
    


    
    return 0
