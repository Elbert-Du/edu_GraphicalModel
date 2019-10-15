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
from privBayes import privBayesSelect


def datavector(df, domain, flatten=True):
    """ return the database in vector-of-counts form """
    if type(domain) is int:
        bins = [range(domain+1)]
    else:
        bins = [range(n+1) for n in domain]
    ans = np.histogramdd(df.values, bins)[0]
    return ans.flatten() if flatten else ans

def transform_data(data, domain, supports):
    #trancate data to have shrinked domain
    df = data.copy()
    newdom = {}
    for col in domain:
        print(col)
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
    #reverse data to have original domain
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

def moments_calibration(sens1, sens2, eps, delta):
    # round1: L2 sensitivity of round1 queries
    # round2: L2 sensitivity of round2 queries
    # works as long as eps >= 0.01; if larger, increase orders
    # if we do not shrink domain, we will not measure round 1 queries, so set round1's sensitivity very small and not compute its sigma
    orders = range(2, 4096)

    def obj(sigma):
        rdp1 = compute_rdp(1.0, sigma/sens1, 1, orders)
        rdp2 = compute_rdp(1.0, sigma/sens2, 1, orders)
        rdp = rdp1+rdp2  
        if(sens1<=0.000000001):
            rdp = rdp2 
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

def prior_err(variance, q_list, domain):
    #conpute confidence interval of selected queries
    m_list=[]
    for q in q_list:
        m=np.prod([domain[a] for a in q])
        m_list.append(m)
    def find(value, m):
        c=np.sqrt(2/np.pi)*value*np.exp(-m**2/(2*(value**2)))
        return c-0.05
    workload_confidence={}
    for i,m in enumerate(m_list):
        value=optimize.bisect(find, 1, m, args=(m))
        intervals=m*np.sqrt(variance)/value
        workload_confidence[str(q_list[i])]=intervals
    return workload_confidence
    

def r_to_python(queries_r, columns_names):
    #only for new query selection method. just to change the format of query list
    selected_queries=[]
    for marg in queries_r:
        temp=tuple(marg.split(','))
        marg_names=[]
        for t in temp:
            marg_names.append(columns_names[int(t)-1])
        selected_queries.append(tuple(marg_names))
    return selected_queries

class Match3(Mechanism):
    #Please set iters as 10000 
    def __init__(self, dataset, specs, domain, mapping, delta,save="out.csv",iters=10000, warmup=False, is_encoded = True):
        Mechanism.__init__(self, dataset, specs, domain, mapping)
        self.iters = iters
        self.warmup = warmup
        self.elimination_order = None
        self.save=save
        self.is_encoded = is_encoded
        self.data=self.load_data(is_encoded=self.is_encoded)
        self.delta=delta
        self.measurements = []
        self.supports = {}
        for col in self.data.columns:
            self.supports[col] = np.asarray([True for j in range(self.domain[col])])
        
    def shrink_domain(self,epsilon, delta=2.2820610e-12, bound = 0):
        # measure all one-way marginals and shrink domain according to noisy measurements
        data=self.data
        self.round1=list(self.column_order)
        self.delta=delta
        sigma = moments_calibration(1.0, 1.0, epsilon, delta)
        self.sigma = sigma
        print('NOISE LEVEL:', sigma)
        supports = {}
        for i,col in enumerate(self.round1):
            if self.domain[col] <= bound:
                supports[col] = np.asarray([True for j in range(self.domain[col])])
                del(self.round1[i])
        # round1 measurements will be weighted to have L2 sensitivity 1       
        weights = np.ones(len(self.round1))
        weights /= np.linalg.norm(weights) # now has L2 norm = 1
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
            #If the column has fewer possible values than the bound we just leave it be.  
            supports[col] = sup       
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
            
        self.supports=supports
        data,new_domain=transform_data(data,self.domain,supports)
        self.domain=new_domain
        self.data=data
        return data, new_domain

            
    def measure(self,round2, from_r=False):
        # measure selected queries, which are round2 queries
        print("selected queries:",round2)
        if from_r:
            round2 = r_to_python(round2, list(self.column_order))
        self.round2 = round2  #round2 is a query list[]
        # round2 measurements will be weighted to have L2 sensitivity 1
        # perform round 2 measurments over compressed domain
        weights = np.ones(len(self.round2))
        weights /= np.linalg.norm(weights) # now has L2 norm = 1
        for proj, wgt in zip(self.round2, weights):
            #########################
            ## Noise-addition step ##
            #########################
            indices = itemgetter(*proj)(self.domain)
            hist = datavector(self.data[list(proj)], indices)
            Q = matrix.Identity(hist.size)
            noise = self.sigma*np.random.randn(Q.shape[0])
            y = wgt*Q.dot(hist) + noise
            self.measurements.append( (Q, y/wgt, 1.0/wgt, proj) )

    def postprocess(self):
        #use noisy measurements to fit PGM inference
        #and generate synthetic data
        iters = self.iters
        domain = self.domain
        temp_domain = Domain.fromdict(domain)
        engine = FactoredInference(temp_domain,
                                    structural_zeros=None,
                                    iters=iters,
                                    log=True,
                                    warm_start=False,
                                    elim_order=self.elimination_order)
        self.engine = engine
        engine.estimate(self.measurements)

        self.synthetic = self.engine.model.synthetic_data()
        self.synthetic = reverse_data(self.synthetic, self.supports)
#        print("postprocess:",self.synthetic.df)


    def privbayes_query_selection(self,eps,theta,seed):
        #cite: changed from privbayes methods in Ektelo
        #use privbayes to select a set of marginal queries
        domain=self.domain
        config = ''
        for a in list(self.data.columns):
            values = [str(i) for i in range(domain[a])]
            config += 'D ' + ' '.join(values) + ' \n'
        config = config.encode('utf-8')
        values = np.ascontiguousarray(self.data.values.astype(np.int32))     
        ans = privBayesSelect.py_get_model(values, config, eps, theta, seed)
        ans = ans.decode('utf-8')[:-1]
        projections = []
        for m in ans.split('\n'):
            p = [list(self.data.columns)[int(a)] for a in m.split(',')[::2]]
            projections.append(tuple(p))
        return projections
    
def default_params():
    """
    Return default parameters to run this program

    :returns: a dictionary of default parameter settings for each command line argument
    """
    params = {}
    params['dataset'] = pd.read_csv('competitor_pack/data/fire-data-2.csv', nrows = 10000)
    
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
    iters = 750 #7500
    
    mech = Match3(args.dataset, args.specs, args.domain, args.mapping, args.delta, args.save, iters=iters, warmup=True)
    mech.shrink_domain(args.epsilon/2,args.delta)
    round2 = mech.privbayes_query_selection(eps=args.epsilon/2,seed=0)

    mech.run(round2,from_r =False) #round2 is a query list [   ]
    
