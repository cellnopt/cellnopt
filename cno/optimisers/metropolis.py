import time
import random
from cno.optimisers.core import Results
from cno.optimisers.diagnostics import Diagnostics
import numpy

class Simulator(object):

    def eval_func(self):
        raise NotImplementedError



class MH(Diagnostics):
    """


    ::
    
        from cno import *
        s = Steady(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
        def eval_func(param):
            return s.eval_func(param)
        m = MH()
        m.run(eval_func, len(s.model.reactions), nswap=3)


    """
    def __init__(self, N=1000):
        Diagnostics.__init__(self)
        self.N = N
        self._buffer = {}

    def neighbour_function(self):
        pass

    def _swap(self, x):
        if x == 1:
            return 0
        else:
            return 1

    def swaps(self, bitstring, nswap):
        # indexlist = list that store the random indices that are created 
        # swapbits function swaps the bits of the indices that indexlist has
        assert nswap >0 and nswap<=len(bitstring)
        bit = bitstring[:]  # we do not want to modify the input
        # let us pick up N unique indices to swap bits
        indexlist = random.sample(range(0, len(bit)), nswap)
        # let us do the swaps now
        newbitstring =  self._swapbits(bit, indexlist)
        return newbitstring

    def maxswaps(self, bitstring, nswap):
        """ Returns a modified version of the given bitstring by swapping 1 to N bits 
 
        :param bitstring: the list that corresponds to the model
        :param nswap: the maximum number of bits to swap 
      
          Running an example with nswap=5. The new bitstring can have from 1 to 5 zeros

          >>> bitstring = [1,1,1,1,1,1]
          >>> new_bitstring = swaps(bitstring,5)
          >>> new_bitstring
          [1,1,0,0,1,0]

       """
        new = bitstring[:]
        L = len(bitstring)
        assert n<=L
        n = numpy.random.randint(n+1)
        if n == 0:
            n = 1
        indices = set()
        while len(indices) != nswap:
            r = numpy.random.randint(L)
            indices.add(r)
        for i in indices:
            new[i] = self._swap(bitstring[i])
        return new

    def run(self, eval_func, N, nswap=3, proposal=None):

        self.Nparameters = N

        results = Results(N=N, step=1)

        t1 = time.time()
        self.alpha = []
        if proposal is None:
            proposal_parameter = [1] * self.Nparameters
        else:
            proposal_parameter = proposal[:]
        init_bs = proposal_parameter[:]

        prev_score = eval_func(proposal_parameter)  # compute the score for the initial bitstring
        prev_bs = init_bs[:]
        best_score = prev_score
        results['best_score'] = best_score
        best_parameters = init_bs[:]

        # store the initial values
        results['scores'].append(prev_score)
        results['parameters'].append(prev_bs)

        from easydev import progress_bar
        pb = progress_bar(self.N)

        for i in range(0, self.N):

            proposal_parameter = self.swaps(best_parameters, nswap)

            #tup_param = tuple(proposal_parameter)
            #if tup_param in self._buffer.keys():
            #    proposal_score = self._buffer[tup_param]
            #else:
            #    proposal_score = eval_func(proposal_parameter)
            #    self._buffer[tup_param] = proposal_score
            proposal_score = eval_func(proposal_parameter)

            alpha = prev_score / proposal_score # best score is the smallests one
            # so alpha >1 means new proposal is better

            self.alpha.append(alpha)

            if alpha >=1:
                prev_score = proposal_score
                score = proposal_score
                results['parameters'].append(proposal_parameter)
                prev_bs = proposal_parameter[:]
                accepted = 1
            else:
                r = random.uniform(0,1)
                if r <= alpha:
                    prev_score = proposal_score
                    score = proposal_score
                    # storing results
                    results['parameters'].append(proposal_parameter)
                    prev_bs = proposal_parameter[:]
                    accepted = 1
                else:
                    prev_score = prev_score
                    # storing results
                    score = prev_score
                    results['parameters'].append(prev_bs)
                    accepted = 0
            self.acceptance.append(accepted)
            results['scores'].append(score)

            if score < best_score:
                best_score = score
                best_parameters = proposal_parameter[:]
                results['best_score'] = best_score
            results['best_scores'].append(best_score)
            results['best_score'] = best_score # just for the progres
            pb.animate(i)

        #print best_parameters
        del results['scores'][0] # remove first element to have a length of N value
        del results['parameters'][0] # remove first element to have a length of N value
        results['best_parameters'] = best_parameters[:]
        results['min_index'] = numpy.argmin(results['best_scores']) # store the index of the minimum score from the best_scores list     

        t2 = time.time()
        print "simulation took", t2-t1, "seconds."

        self.results = results.copy()
