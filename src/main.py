from calendar import c
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import sem

class ExamSolver:
    def __init__(self, l=50, nstep=int(1e5), D=1, q=1, p=0.5):
        self.l, self.N = l, l**2
        self.dx, self.D, self.q, self.p = 1, D, q, p
        self.dt = 0.1
        self.nstep = nstep
        # Initialise random number generator. This is the 
        # better way to generate random numbers, with generator PCG64
        self.rng = np.random.default_rng()
        self.a, self.b, self.c = self.get_initial_concentrations()
        # print(f'Array a is \n{self.a}')
        # print(f'Array b is \n{self.b}')
        # print(f'Array c is \n{self.b}')

    def get_initial_concentrations(self):
        a = self.rng.random(size=(self.l, self.l)) / 3
        b = self.rng.random(size=(self.l, self.l)) / 3
        c = self.rng.random(size=(self.l, self.l)) / 3
        return a, b, c
        
    def laplacian_2d(self, arr: np.ndarray):
        """Computes the laplacian for a 2 dimensional array"""
        return (np.roll(arr, 1, axis=0) + np.roll(arr, 1, axis=1) + \
        np.roll(arr, -1, axis=0) + np.roll(arr, -1, axis=1) - 4 * arr) / self.dx**2

    def update_a(self, a, b, c):
        self.a += self.dt * (
            self.D * self.laplacian_2d(a) + self.q * a * (1 - a - b - c) - self.p * a * c
        )

    def update_b(self, a, b, c):
        self.b += self.dt * (
            self.D * self.laplacian_2d(b) + self.q * b * (1 - a - b - c) - self.p * a * b
        )

    def update_c(self, a, b, c):
        self.c += self.dt * (
            self.D * self.laplacian_2d(c) + self.q * c * (1 - a - b - c) - self.p * b * c
        )

    def update_concentrations(self):
        # Take copies of concentration arrays a, b and c so that
        # the updates do not use the newly updated concentrations
        # params = {'a': self.a.copy(), 'b': self.b.copy(), 'c': self.c.copy()}
        current_a, current_b, current_c = self.a.copy(), self.b.copy(), self.c.copy()
        self.update_a(current_a, current_b, current_c)
        self.update_b(current_a, current_b, current_c)
        self.update_c(current_a, current_b, current_c)

    def get_tau_value(self, a_ij, b_ij, c_ij):
        d_ij = 1 - a_ij - b_ij - c_ij
        values = [d_ij, a_ij, b_ij, c_ij]
        # print(values)
        # max_val = max(values)
        idx = values.index(max(values))
        # print(f'Max values is in elelmnet {idx}')
        # print(idx)
        return idx

    def get_field_tau(self):
        tau = np.zeros(shape=(self.l, self.l))
        for row in range(self.l):
            for col in range(self.l):
                tau[row, col] = self.get_tau_value(
                    self.a[row, col],
                    self.b[row, col],
                    self.c[row, col]
                    )
        return tau


    def visualise(self):
        tau = self.get_field_tau()
        plt.imshow(tau, cmap='plasma', vmin=0, vmax=3)
        plt.colorbar()

        for epoch in range(self.nstep):
            # Update concentrations according to the PDEs
            self.update_concentrations()
            tau = self.get_field_tau()

            if epoch % 10 == 0:
                print(f"Epoch number {epoch}")
                plt.cla()
                plt.imshow(tau, animated=True, cmap='plasma', vmin=0, vmax=3)
                plt.title(r'$\tau$ field at time ' + str(epoch))
                plt.draw()
                plt.pause(0.0001)

    def get_fraction(self, tau, concentration_type: int):
        return np.sum(tau == concentration_type) / self.N

    def get_fraction_of_concentrations(self):
        times, frac_a, frac_b, frac_c = 0, [], [], []
        tau = self.get_field_tau()
        while not self.is_absorbing_state(tau):
            self.update_concentrations()
            tau = self.get_field_tau()
            frac_a.append(self.get_fraction(tau, concentration_type=1))
            frac_b.append(self.get_fraction(tau, concentration_type=2))
            frac_c.append(self.get_fraction(tau, concentration_type=3))
            times += 1

        data = {'time': list(range(1, times+1)), 'fraction_a': frac_a, 'fraction_b': frac_b, 'fraction_c': frac_c}
        pd.DataFrame.from_dict(data).to_csv('results/fraction_of_points_question_b.csv', index=False)


    def is_absorbing_state(self, tau):
        return np.all(tau == 1) or np.all(tau == 2) or np.all(tau ==3)

    def get_time_to_absorption(self):
        n = int(int(1e3) / self.dt)
        for time in range(n):
            self.update_concentrations()
            tau = self.get_field_tau()
            if self.is_absorbing_state(tau):
                # this is time to absorption
                return time
        # If after n simulations there is no absorbing state
        # disregard simulation
        return 0

    def get_average_absorption_time(self, num_simulations=10):
        absorption_times = []
        while len(absorption_times) < num_simulations:
            # reset new simulation:
            self.a, self.b, self.c = self.get_initial_concentrations()
            # obtain absorption time for current simulation
            absorption_time = self.get_time_to_absorption()
            # if absorption_time is zero disregard_simulations
            # otherwise add to our list
            if absorption_time != 0:
                absorption_times.append(absorption_time)

        # save absorption time for simulations:
        data = {'simulation': list(range(1, num_simulations+1)), 'absorption_time': absorption_times}
        pd.DataFrame.from_dict(data).to_csv('results/simulations_and_absorption_times.csv', index=False)

        print(f'Average absorption time and standard error on the mean (format mu ± std. err.) is:\n\
            {np.mean(absorption_times)} ± {sem(absorption_times)}')

        print(f'Standard error on the mean: {np.std(absorption_times)/np.sqrt(len(absorption_times))}')

    def get_value_of_a_field_in_two_points(self):
        # get randomly the two points in the grid from which
        # we are going to record their values
        row1, row2, col1, col2 = self.rng.integers(0, self.l, size=4)
        point1_vals, point2_vals = [], []
        for time in range(5000):
            value1 = self.a[row1, col1]
            value2 = self.a[row2, col2]
            point1_vals.append(value1)
            point2_vals.append(value2)

            # Update concentration field
            self.update_concentrations()
            # no need to update tau-field since we do not visualise
            # anything for this question

        data = {'time': list(range(5000)), 'value1': point1_vals, 'value2': point2_vals}
        pd.DataFrame.from_dict(data).to_csv('results/two_points_data_question_e.csv', index=False)
        
    def question_f_single_row_for_a_given_time(self, tau, row):
        xs = defaultdict(int)
        # take a single row
        tau_row = tau[row]
        for idx in range(self.l):
            # we are looping from the 1st element to the last element
            curr_el = tau_row[idx]
            for distance in range(1, 26):
                # here we loop from the right element of idx1 to the following 25.
                # Remember periodic bc. so when we get to 25 we come back to 1
                if curr_el == tau_row[(idx + distance) % self.l]:
                    # two cells with distance 'distance' have the same field
                    # record value in xs
                    xs[distance] += 1
        return xs

    def question_f_all_rows_for_a_given_time(self, tau):
        # store the values for the comparisons in all rows at a given time
        results = defaultdict(float)
        for row in range(self.l):
            xs = self.question_f_single_row_for_a_given_time(tau, row)
            for distance in range(1, 26):
                # add the results for a single row 
                results[distance] += xs[distance]
        # divide by total number of comparisons to get probabilities
        for distance in range(1, 26):
            results[distance] /= self.N

        return results

    def question_f(self, D=0.5):
        self.D, self.q, self.p = D, 1, 2.5
        average_probabilities_as_distance = defaultdict(float)
        times = 1000
        for time in range(1,times+1):
            print(f'Time at {time}')
            self.update_concentrations()
            tau = self.get_field_tau()
            probs_as_distance = self.question_f_all_rows_for_a_given_time(tau)
            for distance in range(1, 26):
                average_probabilities_as_distance[distance] += probs_as_distance[distance]
            
        # After 1000 timesteps take mean per column
        avg_probabilities = []
        for distance in range(1,26):
            # append average probabilities in order to store them in csv file
            avg_probabilities.append(average_probabilities_as_distance[distance]/times)
        
        data = {'p': avg_probabilities, 'r': list(range(1,26))}
        pd.DataFrame.from_dict(data).to_csv(f'results/questionf.csv_D_{self.D}.csv', index=False)
        

if __name__ == '__main__':
    # from question a to c we use these parameters
    params_until_question_c = {'l': 50, 'D': 1 , 'q': 1, 'p': 0.5}

    # create object of class ExamSolver with the parameters we want
    # solver = ExamSolver(**params_until_question_c)

    # Uncomment this to visualise question a-b-c
    # solver.visualise()

    # This method was used to get the fraction of the concentrations
    # asked in question b
    # solver.get_fraction_of_concentrations()

    # This method was used to get the average absorption time
    # asked in question c
    # solver.get_average_absorption_time()

    # from question d to e we use these parameters
    params_from_question_d = {'l': 50, 'D': 0.5 , 'q': 1, 'p': 2.5}
    # This was used to get the animation for question d
    solver = ExamSolver(**params_from_question_d)
    solver.visualise()

    # This method was used to obtain the value of the a field in two
    # different points for question e
    # solver.get_value_of_a_field_in_two_points()

    # These methods were used to obtain the data required for the plots
    # in question f (probability wrt distance) for the different values
    # of D
    # solver.question_f(0.5)
    # solver.question_f(0.3)
    # solver.question_f(0.4)