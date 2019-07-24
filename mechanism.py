import numpy as np
import pandas as pd
import json
from scipy.stats import norm
import pickle
#from mbi import Dataset, Domain

class Mechanism:
    """ This class is a template for a mechanism with all the boilerplate code
        already implemented.  subclasses should implement three functions:
            setup, measure, and postprocess

        measure is the only function that is allowed to look at the data, and it must
        be privacy-vetted.  All other code should not need to be checked with very much scrutiny
    """
    def __init__(self, dataset, specs, domain, mapping):
        self.mapping = mapping
        self.dataset = dataset
        self.specs = specs
        domain_info = domain
        del(domain_info["columns"])

        # check consistency for codebook information
        #for col in list(domain_info):
        #    if domain_info[col][-1] < self.specs[col]['maxval']:
        #        print('Codebook inconsistent for', col)
        #        del domain_info[col]

        ## look at ground truth data to obtain possible values for state-dependent columns
        #df = pd.read_csv(dataset)
        #for col in ['SEA', 'METAREA', 'COUNTY', 'CITY', 'METAREAD']:
        #    domain_info[col] = sorted(df[col].unique())
        ## done using ground truth data 

        domain = { }
        for col in domain_info:
            domain[col] = len(domain_info[col])

        #domain['INCWAGE_A'] = 52
        #domain['INCWAGE_B'] = 8
        #del domain['INCWAGE']
        #domain['INCWAGE'] = 5002
        #domain['VALUEH'] = 5003
        
        self.domain_info = domain_info 
        self.domain = domain

    def setup(self):
        """ do any setup needed to run the algorithm here """
        pass
    
    def float_map(self, df, col):
        this_min = self.domain_info[col][0]
        this_max = self.domain_info[col][-1]
        step_size = this_max-this_min
        new_df = pd.Series(df[col])
        for i,element in enumerate(df[col]):
            new_df[i] = np.floor((element-this_min)/step_size)
        return(new_df)

    def load_data(self, df = None, is_encoded = False):
        """ load the data and discretize the integer/float attributes """
        #Already discretized in domain.ipynb, but for consistency we map to 0:d-1 for d-1 possible values
        if df is None:
            df = self.dataset
        self.column_order = df.columns
        if not is_encoded:
            for col in self.domain_info:
                if self.specs[col]["type"] == "enum" or self.specs[col]["type"] == "integer":
                    if self.specs[col]["type"] == "integer":
                        df[col].astype('int32', errors = "ignore")
                    vals = self.domain_info[col]
                    mapping = dict(zip(vals, range(len(vals))))
                    df[col] = df[col].map(mapping)

                else:
                    df[col] = self.float_map(df, col)

        #print(df.head())
    

        return df


    
    
    def measure(self):
        """ load the data and measure things about it
        save the measuremnts taken, but do not save the data 
        this is the only function that needs to be vetted for privacy
        """
        pass

    def postprocess(self):
        """ post-process the measurments taken into a synthetic dataset over discrete attributes
        """

    def transform_domain(self, df, mapping):
        """ convert the synthetic discrete data back to the original domain
            and add any missing columns with a default value """
        for col in df:
            #print(mapping[col])
            if col in mapping:
                df[col] = df[col].map(mapping[col])
                
        return df


    def run(self, save=None,round2):
        """ Run the mechanism at the given privacy level and return the synthetic data

        :param epsilon: the privacy budget
        :param delta: privacy parameter
        :param save: location to save the synthetic data
        :return: the synthetic data in the same format as original data
        """
#        self.epsilon = epsilon
#        self.delta = delta
        self.save = save
#        self.setup()
        self.measure(round2)
        self.postprocess()
        self.synthetic.df = self.transform_domain(self.synthetic.df, self.mapping)
        if save is not None:
            self.synthetic.df.to_csv(save, index=False)
        return self.synthetic

if __name__ == '__main__':
    from IPython import embed
    mech = Mechanism()
    df = mech.load_data().df
    mech.synthetic = df
    mech.transform_domain()
    df2 = mech.synthetic
    #embed()
