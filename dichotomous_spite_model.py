"""This code has been written by Pedro Batista Tan, for his major internship of the master program Systems biology and
Bioinformatics, from the Vrije University of Amsterdam (VU) and Universiteit van Amsterdam (UvA). The script has been
adapted from a previous version written by Silvia Espada Burriel. The code has been written with the supervision of
Rutger Hermsen, from the Theoretical biology group of the University of Utrecht."""

"""The code requires setting up parameters for the model simulations. This is done in the first few lines by
 constructing a dictionary with parameter names as keys and lists as values, where every index in the list represents
 a different parameter combination. As currently implemented, the dictionary is built in a for loop with a factorial
 design, using every parameter combination provided in parameter lists.
"""

# IMPORTS ---------------------------------------------------------------------
import sys
import numpy as np
from scipy.stats import norm
from time import time
from datetime import timedelta, datetime
import os
import math
import pathlib

results_folder_name = "Results_uu_dic"
discrete_delta_x = 1  # Used for scaling - how much each grid cell represents of a "physical" space
discrete_delta_t = 1  # Used for scaling - how much each time step represents of a "physical" time

# Create a dictionary of parameters for multiple runs.
# The index of the lists are used to set up the parameter combinations
parameters = {'cost': [],
              'physical_comp_sigma': [],
              'physical_spite_sigma': [],
              'physical_diff': [],
              'spite_expression_prob': [],
              'kd': [],
              'initial_spite_frequency': [],
              'physical_kc': []
              }

# Populate the dictionary with each parameter combination
for physical_diff in [8, 32, 128, "i"]:
    if physical_diff != "i":
        diff_sigma = math.sqrt(2 * physical_diff * discrete_delta_t)/discrete_delta_x
        sweep = [diff_sigma / 4, diff_sigma / 2, diff_sigma, diff_sigma * 2, diff_sigma * 4]
        sweep = list(map(int, sweep))
    else:
        diff_sigma = "i"
        sweep = [2, 4, 8, 16, 32]
    for kd in [0.002]:
        for cost in [0.1]:
            for initial_spite_frequency in [0.5]:
                for spite_expression_prob in [0.001, 0.01, 0.05, 0.2]:
                    for physical_kc in [0.05]:
                        for physical_comp_sigma in sweep:
                            for physical_spite_sigma in sweep:
                                parameters["cost"].append(cost)
                                parameters["physical_comp_sigma"].append(physical_comp_sigma)
                                parameters["physical_spite_sigma"].append(physical_spite_sigma)
                                parameters["physical_diff"].append(physical_diff)
                                parameters["initial_spite_frequency"].append(initial_spite_frequency)
                                parameters["spite_expression_prob"].append(spite_expression_prob)
                                parameters["kd"].append(kd)
                                parameters["physical_kc"].append(physical_kc)


# If directory does not exist, create it. If it exists, check whether results should be overwritten or not.
overwrite_dir_results = False
if not os.path.exists(results_folder_name):
    os.makedirs(results_folder_name)
else:
    if not overwrite_dir_results:
        print("Experiment folder already exists. Please change the name to prevent overwriting or set "
              "overwrite_dir_results to True.")
        sys.exit()

# Perform a simulation and save results for each parameter combination
n_combinations = len(parameters["cost"])
Parameter_sweep_start = time()

for i in range(n_combinations):
    print(f"Parameter combination: {i} / {n_combinations - 1}")
    cumulative_run_time = time() - Parameter_sweep_start
    print(f"Cumulative execution time: {timedelta(seconds=cumulative_run_time)}")

    print("cost", parameters["cost"][i])
    print("physical_comp_sigma", parameters["physical_comp_sigma"][i])
    print("physical_spite_sigma", parameters["physical_spite_sigma"][i])
    print("physical_diff", parameters["physical_diff"][i])
    print("spite_expression_prob", parameters["spite_expression_prob"][i])
    print("initial_spite_frequency", parameters["initial_spite_frequency"][i])
    print("kd", parameters["kd"][i])
    print("physical_kc", parameters["physical_kc"][i])

    # PARAMETERS ------------------------------------------------------------------
    OBS = ""  # Observations saved to log file
    XY = 2 ** 9  # dimensions of the field (2**n)
    N = 13000  # initial number of individuals
    t, t_max = 0, 10_000  # number of iterations - time steps
    past_mean_pheno = 1
    mu = 0  # general mean of the different processes

    initial_spite_frequency = parameters["initial_spite_frequency"][i]  # initial spite carrier frequency

    spite_expression_prob = parameters["spite_expression_prob"][i]  # Chance that carriers burst releasing spite
    scaled_spite_expression_prob = spite_expression_prob * discrete_delta_t

    # probability and step of a mutation event
    pmut = 0.01
    mut_step = 0.5  # large step for the dichotomous model - represents the fixed spite level of carriers

    # growth rate formula parameters
    physical_kc = parameters["physical_kc"][i]  # "maximal" carrying capacity
    kc = physical_kc * discrete_delta_x**2  # scaled carrying capacity
    physical_g0 = 1   # basal growth rate
    g0 = physical_g0 * discrete_delta_t  # scaled basal growth rate
    cost = parameters["cost"][i]  # Spite cost

    # death rate formula parameters
    physical_d0 = 0.1  # basal death rate
    d0 = physical_d0 * discrete_delta_t  # scaled basal death rate
    kd = parameters["kd"][i]  # Tolerance to spite density
    carrier_resistance = True  # determines if spite level is added to the individual's kd

    physical_comp_sigma = parameters["physical_comp_sigma"][i]  # bandwidth of the competition kernel
    comp_sigma = physical_comp_sigma * discrete_delta_x

    physical_spite_sigma = parameters["physical_spite_sigma"][i]   # bandwidth of the spite kernel
    spite_sigma = physical_spite_sigma * discrete_delta_x

    physical_diff = parameters["physical_diff"][i]  # diffusion constant in length^2/time

    if physical_diff == "i":
        diff_sigma = "i"
        scaled_diff_sigma = 0  # In this case this parameter is not used, but is still defined here
    else:
        diff_sigma = math.sqrt(2 * physical_diff * discrete_delta_t)/discrete_delta_x  # standard deviation of diffusion
        scaled_diff_sigma = diff_sigma/math.sqrt(2)  # Scaled SD of diffusion, as there are 2 diffusions per time step
        diff_sigma = int(diff_sigma)  # rounded to an integer (just for the filename)

    fields = max(np.int_((6 * comp_sigma) / XY), 3)  # number of fields to take into account for the kernel estimation.
    # Take 3 if sigma is much lower than XY

    # Determine path to save experiment results. A naming convention with the parameter combinations was used
    folder = pathlib.Path(__file__).parent.absolute()  # get parent - current script directory
    dir_name = f"d{diff_sigma}_c{comp_sigma}_s{spite_sigma}_cost{cost}" \
               f"_kd{kd}_init{initial_spite_frequency}_eprob{spite_expression_prob}" \
               f"_res{carrier_resistance}_kc{physical_kc}"

    path = '%s/%s/%s' % (folder, results_folder_name, dir_name)  # full path

    # If directory does not exist, create it. If it exists, check whether results should be overwritten or not.
    if not os.path.exists(f"{results_folder_name}/{dir_name}"):
        os.makedirs(f"{results_folder_name}/{dir_name}")
    else:
        if not overwrite_dir_results:
            print("Experiment folder already exists within results folder. Please change the experiment name to prevent"
                  "overwriting or set overwrite_dir_results to True.")
            sys.exit()

    # FUNCTIONS -------------------------------------------------------------------


    def first_individuals(N, initial_spite_frequency):
        """Generate the first N individuals at random positions with the initial frequency of carriers"""

        # Coordinates: two columns and N rows of ints with max of XY (not counted)
        x = np.random.randint(XY, size=N)
        y = np.random.randint(XY, size=N)

        # Variables type definition of the structured array
        dt = np.dtype([('x', np.int16), ('y', np.int16), ('p', np.float32), ('w', np.int16), ('d', np.float32)])
        # Pre-definition of the array filled with zeros
        individuals = np.zeros(N, dtype=dt)

        # choose indices to convert to initial carriers
        initial_spite_number = int(N*initial_spite_frequency)
        initial_carrier_indexes = np.random.choice(range(N), size=initial_spite_number, replace=False)

        # substitute the zeros with values
        individuals['x'] = x
        individuals['y'] = y
        individuals['w'] = 1
        # transform individuals into carriers at the chosen indexes, add mut_step to spite ['p'] = 0
        np.add.at(individuals['p'], initial_carrier_indexes, mut_step)

        return individuals


    def report(individuals, t, stats):
        """ At every time-step, save the time, the number of individuals, the spite carrier proportion,
        the standard deviation of the phenotype, the covariance between
        the phenotype and the fitness and the expected value. Save it as a file every
        time-step, re-writing the file, instead of at the end in case it stops prematurely."""

        global past_mean_pheno  # to change it globally and keep track of the last value

        # Calculate the covariance between the spite genotype and fitness (cov_w) or the spite genotype and the
        # spite density (cov_d) through a simplified formula using means
        mean_pheno = np.mean(individuals['p'])
        mean_w = np.mean(individuals['w'])
        mean_pheno_density = np.mean(individuals['d'])
        max_pheno_density = np.max(individuals['d'])
        cov_w = np.mean((individuals['p'] - mean_pheno) * ((individuals['w']) / mean_w - 1))
        cov_d = np.mean((individuals['p'] - mean_pheno) * (individuals['d'] - mean_pheno_density))
        var_p = np.var(individuals['p'])
        spite_carrier_proportion = sum(individuals['p'] > 0)/len(individuals)

        stats['t'][t] = t
        stats['n'][t] = len(individuals)
        stats['mp'][t] = mean_pheno
        stats['pp'][t] = spite_carrier_proportion
        stats['s'][t] = cov_w
        stats['e'][t] = (mean_pheno - past_mean_pheno) - cov_w
        stats['d'][t] = mean_pheno_density
        stats['maxd'][t] = max_pheno_density
        stats['r'][t] = cov_d/var_p

        if t % (t_max / 1000) == 0:
            np.savez('%s/stats' % path, t=stats['t'], n=stats['n'], mp=stats['mp'], pp=stats['pp'],
                     s=stats['s'], e=stats['e'], d=stats['d'], maxd=stats['maxd'], r=stats['r'])

        past_mean_pheno = mean_pheno  # change globally the variable value


    def save_field(file, t):
        """Save, at a specific time, a file with the field of individuals,
        their phenotypes and the density"""
        filename = 'time%i' % t
        print('%s/%s' % (path, filename))
        np.savez('%s/%s' % (path, filename), x=file['x'], y=file['y'], p=file['p'], d=file['d'])


    def transformed_kernel(mu, sigma):
        """Returns (1) the Fourier transform of a matrix containing the weights of the kernel density estimate and
        (2) the kernel weight at position 0,0. A Gaussian kernel is used with mean(mu) and bandwidth(sigma).
        For competition(growth) use the comp_sigma and for spite(death) use the spite_sigma"""

        # start with a vector of zeros with length of the field
        vector = np.zeros(XY)
        for i in range(XY):
            for f in range(-fields, fields + 1):
                # fill each position of the vector with the density of the
                # normal distribution at this position in each field.
                vector[i] += norm.pdf(i + f * XY, loc=mu, scale=sigma)

        # multiply the vector with the transpose of itself to get the kernel matrix
        kernel = np.outer(vector, vector.T)
        # Kernel Normalization
        kernel /= np.sum(kernel)
        # fourier transform of the kernel matrix
        fft_kernel = np.fft.fft2(kernel)

        return fft_kernel, kernel[0, 0]


    def indi_count(individuals):
        """Count the number of individuals at each position. For each position,
        check whether a individual has that coordinates, and if so, add 1"""
        # define the field that will be filled
        field_count = np.int_(np.zeros((XY, XY)))
        x = individuals['x']
        y = individuals['y']
        # Vectorized operation with np "fancy indexing", add 1 to each position for each (x,y) coordinate
        np.add.at(field_count, (x, y), 1)

        return field_count


    def pheno_sum(individuals, burst_carriers):
        """Sum the phenotype at each position where individuals have burst. The index of individuals that burst are used
        to add the value of their phenotypes to the corresponding coordinates."""
        field_pheno = np.zeros((XY, XY))
        spite_x = np.take(individuals['x'], burst_carriers)
        spite_y = np.take(individuals['y'], burst_carriers)
        spite_p = np.take(individuals['p'], burst_carriers)
        # Vectorized operation with "fancy indexing", add corresponding spite to each position for each (x,y) coordinate
        np.add.at(field_pheno, (spite_x, spite_y), spite_p)

        return field_pheno


    def mutation(individuals, nnew):
        """Mutate the phenotype of the new individuals. The mutation happen with
        probability pmut and consist in a alteration of the phenotype with a
        mut_step parameter. For the dichotomous model, spite is set to zero or to the mut_step"""

        prob = np.random.uniform(size=nnew)
        n = 0  # to avoid a nested for loop
        # from the list of individuals, select only for the new ones (end list)
        for i in individuals[(len(individuals) - nnew):]:
            if 0 < prob[n] <= pmut:
                i['p'] = mut_step
            # same probability (equal density) to obtain and lose plasmid
            elif pmut < prob[n] <= 2 * pmut:
                i['p'] = 0
            n += 1

        return individuals


    def growth(individuals):
        """Replicate individuals depending on their growth rate. The growth_rate is calculated based on the competition
        density, the basal growth rate and the metabolic spite cost. Obtain the indices of the lucky individuals to a
        list, copy them and add them to the array. After that, make the new individuals mutate their phenotype."""

        # Obtain the array with the number of individuals at each position and apply the kernel fourier transform.
        count = indi_count(individuals)
        fft_count = np.fft.fft2(count)
        # Multiply them and apply the inverse fourier transform to the result to get the density
        fft_comp_density = fft_kernel_comp * fft_count
        comp_density = (np.fft.ifft2(fft_comp_density)).astype(float)

        # obtain a vector with the competition density for each position, discounting the individual's own contribution
        position_comp_densities = comp_density[individuals['x'], individuals['y']] - position0_comp
        growth_rates = g0 * (1 - position_comp_densities / kc) * (1 - cost * individuals['p'])

        # Get a list of individual indices to replicate
        prob = np.random.uniform(size=len(individuals))
        new = growth_rates > prob
        new = np.nonzero(new)[0]  # the desired array is the first element of a tuple returned from np.nonzero()
        new = new.tolist()

        # Replicate individuals
        individuals = np.append(individuals, individuals[new], axis=0)

        # Add 1 to the individual's fitness that have had offspring
        individuals['w'][new] += 1

        # Apply mutations to newborn individuals
        individuals = mutation(individuals, len(new))
        return individuals


    def death(individuals, t):
        """Remove individuals depending on their death rate and individuals that burst releasing spite. The death_rate
         is calculated with the density of released spite and other parameters. Obtain the indices of the unlucky
         individuals and delete them from the array. A report of all individuals is made before any is removed"""

        # Obtain a list of the carrier indices
        carriers = individuals['p'] > 0  # Value of 1 for carriers, 0 for non carriers
        carriers = np.nonzero(carriers)[0]  # array is the first element of a tuple returned from np.nonzero()
        carriers = carriers.tolist()

        # Determine which carriers should burst and extract a list from the carrier indices
        prob = np.random.uniform(size=len(carriers))
        spite_expression = scaled_spite_expression_prob > prob
        burst_carriers = np.extract(spite_expression, carriers)

        # count phenotype density only of individuals that burst
        if len(burst_carriers) > 0:
            pheno = pheno_sum(individuals, burst_carriers)
        else:
            pheno = np.zeros((XY, XY))

        # Apply the kernel fourier transform, then
        # multiply them and do the inverse transform to the result to get the density
        fft_pheno = np.fft.fft2(pheno)
        fft_pheno_density = fft_kernel_spite * fft_pheno
        pheno_density = (np.fft.ifft2(fft_pheno_density)).astype(float)

        # obtain a vector with the spite density for each position
        position_pheno_density = pheno_density[individuals['x'], individuals['y']]
        individuals['d'] = position_pheno_density

        if kd == 'i':  # infinite resistance to spite
            death_rate = d0
        elif carrier_resistance:
            death_rate = d0 * (1 + position_pheno_density / (kd + individuals['p']))
        else:
            death_rate = d0 * (1 + position_pheno_density / kd)

        # Check if any death rate is above 1
        # warning = death_rate > 1
        # if np.any(warning):
        #     warning = np.nonzero(warning)[0]
        #     warning = warning.tolist()
        #     print(f"Death rate exceeded 1: {death_rate[warning]}")
        #     print(f"Spite Density: {individuals['d'][warning]}")

        # Get a list of individual indices to kill, also including carriers that burst
        prob = np.random.uniform(size=len(individuals))
        idx = death_rate > prob
        if len(burst_carriers) > 0:
            np.add.at(idx, burst_carriers, 1)
        idx = np.nonzero(idx)[0]
        idx = idx.tolist()

        # Subtract 1 from the fitness of individuals that are going to pass away
        # This is done before deleting them in order to also calculate their fitness
        individuals['w'][idx] -= 1

        # Generate the report inside the death() function because we don't want to
        # miss the individuals that are going to die (be deleted)
        report(individuals, t, stats)

        # after the report, delete the individuals that have died.
        individuals = np.delete(individuals, idx, axis=0)

        # Set to 1 all the individual's fitness for the next steps
        individuals['w'] = 1
        return individuals


    def diffusion(individuals):
        if physical_diff == 'i':
            '''Change coordinate (x,y) to random positions'''
            x = np.random.randint(XY, size=len(individuals))
            y = np.random.randint(XY, size=len(individuals))
            individuals['x'] = x
            individuals['y'] = y

            return individuals

        elif physical_diff == 0:
            return individuals

        else:
            '''For each coordinate (x,y) a random step is added from a normal distribution.
            Note that the normal distribution is being discretized by choosing integers from a range based on the 
            diffusion probabilities. After adding each step, the boundaries are redefined (%) in order to make them
            periodic, so if X>XY: X becomes the remainder of X/XY, if X<XY: X=XY'''

            delta_x = np.random.choice(diff_range, size=len(individuals), p=diff_prob)
            delta_y = np.random.choice(diff_range, size=len(individuals), p=diff_prob)

            individuals['x'] += delta_x
            individuals['y'] += delta_y
            # periodic field boundaries
            individuals['x'] %= XY
            individuals['y'] %= XY

            return individuals


    def save_log_file(path):
        # Save parameters used in log file
        run_start_time = datetime.now().strftime("%d-%m-%Y %Hh%M")
        log_save_name = f"parameters_log.txt"
        path = f"{path}/{log_save_name}"
        file = open(path, "w")
        file.write(f"Script name: {os.path.basename(__file__)}\n")
        file.write(f"run start GMT: {run_start_time}\nXY: {XY}\nN: {N}\nt_max: {t_max}\nmu: {mu}\n")
        file.write(f"discrete_delta_x: {discrete_delta_x}\ndiscrete_delta_t: {discrete_delta_t}\n")
        file.write(f"physical_diff: {physical_diff}\nphysical_comp_sigma: {physical_comp_sigma}\n"
                   f"physical_spite_sigma: {physical_spite_sigma}\n")
        file.write(f"diff_sigma: {diff_sigma}\nscaled_diff_sigma: {scaled_diff_sigma}\ncomp_sigma: {comp_sigma}\n")
        file.write(f"spite_sigma: {spite_sigma}\ninitial_spite_frequency: {initial_spite_frequency}\n"
                   f"pmut: {pmut}\nmut_step: {mut_step}\nfields: {fields}\n")
        file.write(f"physical_g0: {physical_g0}\ng0: {g0}\nphysical_kc: {physical_kc}\nkc: {kc}\n"
                   f"cost: {cost}\nphysical_d0: {physical_d0}\nd0: {d0}\nkd: {kd}\n"
                   f"spite_expression_prob: {spite_expression_prob}\ncarrier_resistance: {carrier_resistance}\n")
        file.write(f"OBS: {OBS}\n")
        file.close()


    # BEFORE START ----------------------------------------------------------------
    run_start = time()

    # Calculate only once the kernels
    fft_kernel_comp, position0_comp = transformed_kernel(mu, comp_sigma)
    fft_kernel_spite, position0_spite = transformed_kernel(mu, spite_sigma)

    # Set range and probabilities for integer displacements in diffusion steps
    diff_bounds = int(scaled_diff_sigma * 6)
    diff_range = np.arange(-diff_bounds, diff_bounds)
    xU, xL = diff_range + 0.5, diff_range - 0.5
    diff_prob = norm.cdf(xU, scale=scaled_diff_sigma) - norm.cdf(xL, scale=scaled_diff_sigma)
    diff_prob = diff_prob / diff_prob.sum()  # normalize the probabilities so their sum is 1

    # Generate array for statistics
    dt_s = np.dtype([('t', np.int32), ('n', np.int32), ('mp', np.float32), ('pp', np.float32),
                     ('s', np.float32), ('e', np.float32), ('d', np.float32), ('maxd', np.float32),
                     ('r', np.float32)])  # Variables type definition
    stats = np.zeros((t_max + 1, 1), dtype=dt_s)

    # Generate the first individuals
    individuals = first_individuals(N, initial_spite_frequency)

    save_log_file(path)

    # CODE ------------------------------------------------------------------------
    # Iterate t_max times and make the individuals die, grow and diffuse.
    # If there are no more individuals, the loop stops.

    while t <= t_max and len(individuals) != 0:
        if t % (t_max/10) == 0:
            save_field(individuals, t)
        individuals = death(individuals, t)
        individuals = diffusion(individuals)
        individuals = growth(individuals)
        individuals = diffusion(individuals)

        if t % (t_max/10) == 0:
            print(f"time: {t}, spite_carrier_proportion: {sum(individuals['p'] > 0)/len(individuals)},"
                  f"n_individuals: {len(individuals)}")
        t += 1

    run_end = time()
    run_time = run_end - run_start
    print(f"Total execution time: {timedelta(seconds=run_time)}\n")
