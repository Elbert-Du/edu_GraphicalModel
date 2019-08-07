from mbi import Dataset
import numpy as np
import benchmarks
#data = Dataset.load('../data/adult.csv', '../data/adult-domain.json')

def err(true, est):
    return np.sum(np.abs(true - est)) / true.sum()


data, workload = benchmarks.adult_benchmark()
total=data.df.shape[0]
domain="/home/home5/yz488/private-pgm/data/adult-domain.json"
#privbayesdata='/home/home5/yz488/private-pgm/examples/0.5privbayes_'
#dualquerydata='/home/home5/yz488/private-pgm/examples/dualquery_'
gmdata='/home/home5/yz488/edu_GraphicalModel/0.25pg+pb '


#pbee=[]
#dqee=[]
gmee=[]
for i in range(10):
    '''
    pb_path=privbayesdata
    pb_path+=str(i)
    pb_path+=".csv"
    print(pb_path)
    syn_data_privbayes = Dataset.load(pb_path, domain)

    dq_path=dualquerydata
    dq_path+=str(i)
    dq_path+=".csv"
    print(dq_path)
    syn_data_dualquery= Dataset.load(dq_path, domain)
    '''
    gm_path=gmdata
    gm_path+=str(i+1)
    gm_path+=" .csv"
    print(gm_path)
    syn_data_r= Dataset.load(gm_path, domain)
    

   # err_pb = []
   # err_dq = []
    err_r = []
    print("ss")
    for p, W in workload:
        true = W.dot(data.project(p).datavector())
        #    print(data.project(p).datavector())
    #    pb = W.dot(syn_data_privbayes.project(p).datavector())
        #   print(syn_data_privbayes.project(p).datavector())
    #    dq_data=syn_data_dualquery.project(p).datavector()
    #    dq_data*=total/dq_data.sum()
     #   dq = W.dot(dq_data)
        #  print(syn_data_dualquery.project(p).datavector())
        r = W.dot(syn_data_r.project(p).datavector())
        # print(syn_data_r.project(p).datavector())
      #  err_pb.append(err(true, pb))
      #  err_dq.append(err(true, dq))
        err_r.append(err(true, r))
   # pbee.append( np.mean(err_pb))
   # dqee.append(np.mean(err_dq))
    gmee.append(np.mean(err_r))


#print('Error of PrivBayes    : %.3f' % np.mean(err_pb))
#print('Error of DualQuery: %.3f' % np.mean(err_dq))
#print('Error of r    : %.3f' % np.mean(err_r))
#print(pbee)
#print(np.mean(pbee))
#print(dqee)
#print(np.mean(dqee))
print(gmee)
print(np.mean(gmee))
