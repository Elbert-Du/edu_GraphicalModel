from mechanism import Mechanism
import mbi
from mbi import Domain, Factor, FactoredInference, GraphicalModel, Dataset
import matrix
import argparse
import numpy as np
import pandas as pd
from scipy import sparse, optimize
from privacy.analysis.rdp_accountant import compute_rdp, get_privacy_spent
from functools import reduce
import os
import json
from operator import itemgetter

def datavector(df, domain, flatten=True):
    """ return the database in vector-of-counts form """
    if type(domain) is int:
        bins = [range(domain+1)]
    else:
        bins = [range(n+1) for n in domain]
    ans = np.histogramdd(df.values, bins)[0]
    return ans.flatten() if flatten else ans

def transform_data(data, domain, supports):
    df = data.copy()
    newdom = {}
    for col in domain:
        support = supports[col]
        size = support.sum()
        newdom[col] = int(size)
        if size < support.size:
            newdom[col] += 1
        mapping = {}
        idx = 0
        for i in range(support.size):
            mapping[i] = size
            if support[i]:
                mapping[i] = idx
                idx += 1
        assert idx == size
        df[col] = df[col].map(mapping)
    return (df, newdom)

def reverse_data(data, supports):
    df = data.df.copy()
    newdom = {}
    for col in data.domain:
        support = supports[col]
        mx = support.sum()
        newdom[col] = int(support.size)
        idx, extra = np.where(support)[0], np.where(~support)[0]
        mask = df[col] == mx
        if extra.size == 0:
            pass
        else:
            df.loc[mask, col] = np.random.choice(extra, mask.sum())
        df.loc[~mask, col] = idx[df.loc[~mask, col]]
    newdom = Domain.fromdict(newdom)
    return Dataset(df, newdom)

def moments_calibration(round1, round2, eps, delta):
    # round1: L2 sensitivity of round1 queries
    # round2: L2 sensitivity of round2 queries
    # works as long as eps >= 0.01; if larger, increase orders
    orders = range(2, 4096)

    def obj(sigma):
        rdp1 = compute_rdp(1.0, sigma/round1, 1, orders)
        rdp2 = compute_rdp(1.0, sigma/round2, 1, orders)
        rdp = rdp1 + rdp2
        privacy = get_privacy_spent(orders, rdp, target_delta=delta)
        return privacy[0] - eps + 1e-8
    low = 1.0
    high = 1.0
    while obj(low) < 0:
        low /= 2.0
    while obj(high) > 0:
        high *= 2.0
    sigma = optimize.bisect(obj, low, high)
    assert obj(sigma) - 1e-8 <= 0, 'not differentially private' # true eps <= requested eps
    return sigma

class Match3(Mechanism):

    def __init__(self, dataset, specs, domain, mapping, iters=1000, weight3=1.0, warmup=False):
        #domain = json.load(open("domain.json"))
        Mechanism.__init__(self, dataset, specs, domain, mapping)
        self.iters = iters
        self.weight3 = weight3
        self.warmup = warmup
        self.elimination_order = None
        self.mapping = mapping

        #Qa = np.zeros((5,52))
        #Qa[0,0] = 1
        #Qa[1,1:10] = 1
        #Qa[2,10:20] = 1
        #Qa[3,20:51] = 1
        #Qa[4,51] = 1

        #self.Q_INCWAGE = Qa
    
    def setup(self):
        self.round1 = list(self.domain)
        selected_queries = []
        a = open("queries.txt", 'r')
        for line in a.readlines():
            temp = tuple(line.strip('\n').split('_'))
            selected_queries.append(temp)
        self.round2 = selected_queries
        print(selected_queries)
                                        
    def measure(self):
        data = self.load_data(is_encoded = True)
        #print(data.head())
        # round1 and round2 measurements will be weighted to have L2 sensitivity 1
        sigma = moments_calibration(1.0, 1.0, self.epsilon, self.delta)
        print('NOISE LEVEL:', sigma)

        weights = np.ones(len(self.round1))
        #weights[self.round1.index('INCWAGE_A')] *= 2.0
        weights /= np.linalg.norm(weights) # now has L2 norm = 1

        supports = {}
  
        self.measurements = []
        for col, wgt in zip(self.round1, weights):
            ##########################
            ### Noise-addition step ##
            ##########################
            proj = (col,)
            hist = np.asarray(data[col].value_counts())
            print(hist)
            noise = sigma*np.random.randn(hist.size)
            y = wgt*hist + noise
          
            #####################
            ## Post-processing ##
            #####################

            sup = y >= 3*sigma

            supports[col] = sup
            print(col, self.domain[col], sup.sum())
            #print(col, len(self.domain[col]), sup.sum())

            if sup.sum() == y.size:
                y2 = y
                I2 = matrix.Identity(y.size)
            else:
                y2 = np.append(y[sup], y[~sup].sum())
                I2 = np.ones(y2.size)
                I2[-1] = 1.0 / np.sqrt(y.size - y2.size + 1.0)
                y2[-1] /= np.sqrt(y.size - y2.size + 1.0)
                I2 = sparse.diags(I2)

            self.measurements.append( (I2, y2/wgt, 1.0/wgt, proj) )

        self.supports = supports
        # perform round 2 measurments over compressed domain
        data, new_domain = transform_data(data,self.domain, supports)
        self.domain = new_domain
        temp_list = []
        for cl in self.round2:
            temp = 1
            for col in cl:
                temp = temp*self.domain[col]
            if temp < 1e6:
                temp_list.append(cl)
        self.round2 = temp_list
                
        #self.round2 = [cl for cl in self.round2 if len(self.domain[cl] < 1e6]
        #weights = np.ones(len(self.round2))
        #weights[self.round2.index(('SEX','CITY','INCWAGE_A'))] *= self.weight3
        #weights[self.round2.index(('SEX','CITY'))] *= 2.0
        #weights[self.round2.index(('SEX','INCWAGE_A'))] *= 2.0
        #weights[self.round2.index(('CITY','INCWAGE_A'))] *= 2.0
        #weights /= np.linalg.norm(weights) # now has L2 norm = 1
   
        for proj, wgt in zip(self.round2, weights):
            #########################
            ## Noise-addition step ##
            #########################
            indices = itemgetter(*proj)(self.domain)
            print(proj, indices)
            hist = datavector(data[list(proj)], indices)
            Q = matrix.Identity(hist.size)

            noise = sigma*np.random.randn(Q.shape[0])
            y = wgt*Q.dot(hist) + noise
            self.measurements.append( (Q, y/wgt, 1.0/wgt, proj) )

    def postprocess(self):
        iters = self.iters
        domain = self.domain
        temp_domain = Domain.fromdict(domain)
        engine = FactoredInference(temp_domain,
                                    structural_zeros=None,
                                    iters=500,
                                    log=True,
                                    warm_start=True,
                                    elim_order=self.elimination_order)
        self.engine = engine
        cb = mbi.callbacks.Logger(engine)

        if self.warmup:
            engine._setup(self.measurements, None)
            oneway = {}
            for i in range(len(self.round1)):
                p = self.round1[i]
                y = self.measurements[i][1]
                y = np.maximum(y, 1)
                y /= y.sum()
                print(temp_domain.project(p), p, y.shape)
                oneway[p] = Factor(temp_domain.project(p), y)
            marginals = {}
            for cl in engine.model.cliques:
                marginals[cl] = reduce(lambda x,y: x*y, [oneway[p] for p in cl])

            theta = engine.model.mle(marginals)
            engine.potentials = theta
            engine.marginals = engine.model.belief_propagation(theta)

        checkpt = self.save[:-4] + '-checkpt.csv'
        for i in range(self.iters // 500):
            
            engine.infer(self.measurements, engine='MD', callback=cb)

            if i % 4 == 3:
                self.synthetic = engine.model.synthetic_data()
                self.synthetic = reverse_data(self.synthetic, self.supports)
                self.synthetic_df = self.transform_domain(self.synthetic.df, self.mapping)
                self.synthetic_df.to_csv(checkpt, index=False)
   
        if os.path.exists(checkpt):
            os.remove(checkpt)

        self.synthetic = engine.model.synthetic_data()
        self.synthetic = reverse_data(self.synthetic, self.supports)

def default_params():
    """
    Return default parameters to run this program

    :returns: a dictionary of default parameter settings for each command line argument
    """
    params = {}
    params['dataset'] = pd.read_csv('competitor_pack/data/fire-data.csv')
    
    params['specs'] = json.load(open('competitor_pack/data/fire-data-specs.json'))
    params['epsilon'] = 1.0
    params['delta'] = 2.2820544e-12
    params['save'] = 'out.csv'
    params['domain'] = json.load(open('domain.json'))
    mapping = json.load(open('competitor_pack/data/fire-data-specs-mapping.json'))
    for col in mapping:
        for key in list(mapping[col]):
            mapping[col][int(key)] = mapping[col][key]
            del mapping[col][key]
    params['mapping'] = mapping
    return params

if __name__ == '__main__':

    description = ''
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(description=description, formatter_class=formatter)
    parser.add_argument('--dataset', help='dataset csv file')
    parser.add_argument('--specs', help='specs json file')
    parser.add_argument('--epsilon', type=float, help='privacy parameter')
    parser.add_argument('--delta', type=float, help='privacy parameter')
    parser.add_argument('--domain', help='domain file')
    parser.add_argument('--save', help='path to save synthetic data to')
    parser.add_argument('--mapping', help='mapping back to original values')

    parser.set_defaults(**default_params())
    args = parser.parse_args()

    if args.epsilon <= 0.3:
        iters = 7500
        weight3 = 8.0
    elif args.epsilon >= 4.0:
        iters = 10000
        weight3 = 4.0
    else:
        iters = 7500
        weight3 = 6.0

    mech = Match3(args.dataset, args.specs, args.domain, args.mapping, iters=iters, weight3=weight3, warmup=True)

    mech.run(args.epsilon, args.delta, args.save)