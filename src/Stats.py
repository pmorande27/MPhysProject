import numpy as np
import matplotlib.pyplot as plt
class Stats(object):

    def __init__(self, measurements) -> None:

        self.measurements =np.array( measurements)
        

    
    def autocorrelation(self, cutoff, plot= False):

        steps = [i for i in range(cutoff)]

        measurements = self.measurements

        average = np.average(measurements)

        sigma_sq = np.var(measurements,ddof=1)

        #sigma_sq = (average_of_sq-average**2)

        results = [0 for i in range(cutoff)]

        for step in steps:

            result = 0

            for i in range(len(self.measurements)-step):

                result += (measurements[i]-average)*(measurements[i+step]-average)
            
            result = (result/(len(self.measurements)-step))/ sigma_sq
            
            results[step] = result 
        
        if plot:

            plt.plot(results,'o')

            plt.show()
        
        return results
        
    def integrated_autoccorelation(self, cutoff):
        a = self.autocorrelation(cutoff=cutoff)[1:]

        return 2*(1/2 + sum(self.autocorrelation(cutoff=cutoff)[1:]))
    

    def estimate(self, start=0):

        measurements = self.measurements
        average = np.average(measurements)
        cutoff = start

        average_of_sq = np.average([m**2 for m in measurements])

        sigma_sq = np.var(measurements,ddof=1)
        integrated_autoccorelation, error_IAT = Stats.autocorrelator(self.measurements)
        #print(integrated_autoccorelation)
        true_variance = sigma_sq*integrated_autoccorelation
        #print(integrated_autoccorelation)

        error = np.sqrt(true_variance/len(measurements))

        """while True:



            integrated_autoccorelation = self.integrated_autoccorelation(cutoff=cutoff)

            true_variance = sigma_sq*integrated_autoccorelation

            error = np.sqrt(true_variance/len(measurements))
            print(cutoff,integrated_autoccorelation)

            if cutoff >= 4*integrated_autoccorelation:
                break
            cutoff+= 1"""
    


        return average, error
    @staticmethod
    def get_measurements(file_name, observable_name, observable):
        lattices = np.load(file_name)
        name = file_name.split('.')[0] + " " + observable_name
        results = [observable(config) for config in lattices]
        np.savetxt(name,results)
    @staticmethod
    def auto_window(IATs, c):
        '''Windowing procedure of Caracciolo, Sokal 1986 to truncate the sum for the IAT.

        IATs: array
            integrated autocorrelation time with increasingly late termination of the sum
        c: float
            defines the window width. For correlations that decay exponentially, a value of 4 of 5 is conventional

        Returns index for array IATs, representing the IAT estimate from the windowing procedure
        ''' 
        ts = np.arange(len(IATs)) # all possible separation endpoints
        m =  ts < c * IATs # first occurrence where this is false gives IAT 
        if np.any(m):
            return np.argmin(m)
        return len(IATs) - 1
    @staticmethod
    def autocorr_func_1d(x):
        '''Computes the autocorrelation of a 1D array x using FFT and the Wiener Khinchin theorem.
        As FFTs yield circular convolutions and work most efficiently when the number of elements is a power of 2, pad the data with zeros to the next power of 2. 
        '''
        x = np.atleast_1d(x)
        if len(x.shape) != 1:
            raise ValueError("invalid dimensions for 1D autocorrelation function")

        def next_pow_two(n):
            i = 1
            while i < n:
                i = i << 1 # shifts bits to the left by one position i.e. multiplies by 2 
            return i
        
        n = next_pow_two(len(x))

        # Compute the FFT and then the auto-correlation function
        f = np.fft.fft(x - np.mean(x), n=2 * n)
        acf = np.fft.ifft(f * np.conjugate(f))[: len(x)].real
        acf /= 4 * n
        # normalize to get autocorrelation rather than autocovariance
        acf /= acf[0]

        return acf
    @staticmethod
    def autocorr_new(y, c=4.0):
        f = np.zeros(y.shape[1])
        for yy in y:
            f += Stats.autocorr_func_1d(yy)
        f /= len(y)
        taus = 2.0 * np.cumsum(f) - 1.0
        window = Stats.auto_window(taus, c)
        
        return np.mean(y),np.sqrt(taus[window]*np.var(y,ddof=1)/len(y))

    @staticmethod
    def autocorrelator_repeats(data, c=4.0):
    
        M, N = data.shape
        ts = np.arange(M)

        # get ACF and its error
        ACFs = np.zeros_like(data)
        for i in range(N):
            ACFs[:,i] = Stats.autocorr_func_1d(data[:,i])

        ACF, ACF_err = np.mean(ACFs, axis=1), np.std(ACFs, axis=1) / np.sqrt(N)

        # get all possible IAT and apply windowing
        IATs = 2.0 * np.cumsum(ACF) - 1.0 # IAT defined as 1 + sum starting from separation=1, but cumsum starts with t=0 for which ACF=1
        break_idx = Stats.auto_window(IATs, c)
        IAT = IATs[break_idx]
        IAT_err = np.sqrt((4*break_idx+2)/data.shape[0]) * IAT #  Madras, Sokal 1988

        return  IAT, IAT_err
    
    @staticmethod
    def autocorrelator(data, c=4.0):
    
        return Stats.autocorrelator_repeats(data.reshape((data.shape[0], 1)))
    @staticmethod
    def correlator(xs, ys):
   
   
        return Stats.correlator_repeats(xs.reshape((xs.shape[0], 1)), ys.reshape((ys.shape[0], 1)))
    
    @staticmethod
    def correlator_repeats(xs, ys):
    
        N, M = xs.shape # length of one data measurement, number of measurements

        # get ACF and its error
        CFs = np.zeros_like(xs)
        for i in range(M):
            CFs[:,i] = Stats.corr_func_1D(xs[:,i], ys[:,i])

        CF, CF_err = np.mean(CFs, axis=1), np.std(CFs, axis=1) / np.sqrt(M)
        
        # correct error by IAT
       
        return CF
    @staticmethod
    def corr_func_1D(x, y):
    
        f = np.fft.fft(x)
        g = np.fft.fft(y)
        cf = np.fft.ifft(f * np.conjugate(g)).real

        return cf