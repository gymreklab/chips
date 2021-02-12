import math
from math import factorial as fac
import copy

class PCR:
    def __init__(self, kept_rate = 0.5):
        self.kept_rate = kept_rate

    def perform_pcr(self, pcr_rounds=10):
        num_samples_init = 1
        prob_init = 1.0
        freq_dict_prev = {num_samples_init: prob_init}

        # run #pcr_rounds of iterations
        for exp_index in range(pcr_rounds):
            freq_dict = {}
            for n_samples in freq_dict_prev:
                if n_samples > 128:
                    continue
                current_prob = freq_dict_prev[n_samples]
                new_samples_freq_list = self._run_one_step(n_samples)
                for (n_new_samples, new_prob) in new_samples_freq_list:
                    n_samples_current = n_samples + n_new_samples
                    if n_samples_current in freq_dict:
                        freq_dict[n_samples_current] += (current_prob * new_prob)
                    else:
                        freq_dict[n_samples_current] = (current_prob * new_prob)
            freq_dict_prev = freq_dict
        return freq_dict


    # private functions
    def _run_one_step(self, n_samples):
        # binomial prob:
        # sum{ nCi * prob^i * (1-prob)^(n-i) } = 1
        freq_list = [] # elem: (n_new_samples, prob)
        for n_new_samples in range(n_samples+1):
            prob = self._nCr(n_samples, n_new_samples) *\
                    math.pow(self.kept_rate, n_new_samples) *\
                    math.pow(1-self.kept_rate, n_samples-n_new_samples)
            freq_list.append((n_new_samples, prob))
        return freq_list

    def _nCr(self, n, r):
        return (fac(n) / (fac(r)*fac(n-r)))

if __name__ == "__main__":
    pcr = PCR(kept_rate = 0.00935/150.0)
    freq_dict = pcr.perform_pcr(pcr_rounds = 150)
    exp = 0
    for n_samples in freq_dict:
        prob = freq_dict[n_samples]
        if n_samples < 10:
            print ("%d\t%s" %(n_samples, prob))
        exp += n_samples * prob
    print ("Expectation: %s" %(exp))
