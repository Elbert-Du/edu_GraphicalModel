import numpy as np
import pandas as pd
import json
from scipy.stats import norm
import pickle

class Mechanism:
    """ This class is a template for a mechanism with all the boilerplate code
        already implemented.  subclasses should implement three functions:
            setup, measure, and postprocess
        measure is the only function that is allowed to look at the data, and it must
        be privacy-vetted.  All other code should not need to be checked with very much scrutiny
    """
    def __init__(self, dataset, specs, domain, mapping):
        for col in mapping:
            for key in list(mapping[col]):
                mapping[col][int(key)] = mapping[col][key]
                del mapping[col][key]
                
        self.mapping= mapping
        self.dataset = dataset
        self.specs = specs
        domain_info = domain

        domain = { }
        for col in domain_info:
            domain[col] = len(domain_info[col])
  
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
        """ load the data 
        discretize the integer/float attributes, if data is not encoded """
        
        if df is None:
            df = self.dataset
            
        for col in list(df):
            if col not in self.mapping:
                del(df[col])
        self.column_order = df.columns
        if not is_encoded:
            for col in self.domain_info:
                if self.specs[col]["type"] == "character":
                    df[col].astype('int32', errors="ignore")
                    vals=self.domain_info[col]
                    mapping=dict(zip(vals, range(len(vals))))
                    df[col] = df[col].map(mapping)

                else:
                    df[col] = self.float_map(df, col)
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
        """ convert the synthetic discrete data back to the original domain in domain file
            and add any missing columns with a default value """
        for col in df:
            if col in mapping:
                df[col] = df[col].map(mapping[col])
                
        return df


    def write_output(self):
        """write synthetic data to output file"""
        # running tests should use the encoded data, so please comment this line of code
        self.synthetic.df = self.transform_domain(self.synthetic.df, self.mapping)
        self.synthetic.df.to_csv(self.save, index=False)
        return self.synthetic
    

'''
    def run(self, round2, from_r =False):
        """ Run the mechanism at the given privacy level and return the synthetic data

        :param epsilon: the privacy budget
        :param delta: privacy parameter
        :param save: location to save the synthetic data
        :return: the synthetic data in the same format as original data
        """
#        self.epsilon = epsilon
#        self.delta = delta
#        self.save = save
#        self.setup()
        self.measure(round2,from_r =False)
        self.postprocess()
        self.synthetic.df = self.transform_domain(self.synthetic.df, self.mapping)
        self.synthetic.df.to_csv(self.save, index=False)
        return self.synthetic
'''
