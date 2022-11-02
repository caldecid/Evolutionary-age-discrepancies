import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(suppress=True, precision=3)
# load BNN package
import sys
sys.path.insert(0, r'/home/silvestr/Documents/npBNN')
import np_bnn as bn
import scipy.stats
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from concurrent.futures import ProcessPoolExecutor

if __name__=="__main__":
    # set random seed
    rseed = 1234
    np.random.seed(rseed)

    f = "ages.training.predictors.txt"
    l = "ages.training.response.txt"

    import pandas as pd
    f_tbl = pd.read_csv(f, delimiter=',')
    l_tbl = pd.read_csv(l, delimiter=',')

    #### FIT MODEL TO LOG10 AGES
    # rescale features
    f_tbl['Estimated_age'] = np.log10(f_tbl['Estimated.age'])
    f_tbl['root_age'] = np.log10(f_tbl['root.age'])
    f_tbl['number_sp'] = np.log10(f_tbl['number.sp'])
    f_tbl['N_sisters'] = np.log10(f_tbl['N_sisters'])
    # rescale labels
    l_tbl = np.log10(l_tbl)

    seed=rseed

    dat = bn.get_data(f_tbl,
                      l_tbl,
                      seed=seed,
                      testsize=0.1, # 10% test set
                      all_class_in_testset=0,
                      cv=0, # cross validation (1st batch; set to 1,2,... to run on subsequent batches)
                      header=True, # input data has a header
                      from_file=False,
                      instance_id=0,
                      randomize_order=True,
                      label_mode="regression")

    raw_mse = np.mean((dat['test_labels'].flatten() - dat['test_data'][:,0].flatten())**2)

    # set up the BNN model
    bnn_model = bn.npBNN(dat,
                         n_nodes = [8,2],
                         estimation_mode="regression",
                         actFun = bn.ActFun(fun="tanh"),
                         p_scale=1,
                         use_bias_node=3)

    # set up the MCMC environment
    mcmc = bn.MCMC(bnn_model,
                   update_ws=[0.025,0.025, 0.05],
                   update_f=[0.005,0.005,0.05],
                   n_iteration=10000000,
                   sampling_f=1000,
                   print_f=10000,
                   n_post_samples=100,
                   adapt_f=0.3,
                   estimate_error=True)

    # initialize output files
    logger = bn.postLogger(bnn_model, filename="bnn%s" % seed, log_all_weights=0)


    #### MCMCMC
    # set up the BNN model

    n_chains = 10
    rseeds = np.random.choice(range(1000,9999), n_chains, replace=False)


    bnnList = [bn.npBNN(dat,
                        n_nodes = [8,2],
                        estimation_mode="regression",
                        actFun = bn.ActFun(fun="tanh"),
                        p_scale=1,
                        use_bias_node=3,
                        seed=rseeds[i], init_std=0.1)
               for i in range(n_chains)]


    temperatures = np.linspace(0.8, 1, n_chains)
    mcmcList = [bn.MCMC(bnnList[i],
                        temperature=temperatures[i],
                        mcmc_id=i, randomize_seed=True,
                        update_ws=[0.025,0.025, 0.05],
                        update_f=[0.005,0.005,0.05],
                        n_iteration=100,
                        sampling_f=100,
                        print_f=100,
                        n_post_samples=100,
                        adapt_f=0.3,
                        estimate_error=True)
                for i in range(n_chains)]


    singleChainArgs = [[bnnList[i],mcmcList[i]] for i in range(n_chains)]
    n_iterations = 100
    # initialize output files
    logger = bn.postLogger(bnnList[0], filename="BNNMC3", log_all_weights=0)

    def run_single_mcmc(arg_list):
        [bnn_obj, mcmc_obj] = arg_list
        for i in range(n_iterations-1):
            mcmc_obj.mh_step(bnn_obj)
        bnn_obj_new, mcmc_obj_new = mcmc_obj.mh_step(bnn_obj, return_bnn=True)
        return [bnn_obj_new, mcmc_obj_new]

    for mc3_it in range(500):
        with ProcessPoolExecutor(max_workers=n_chains) as pool:
            singleChainArgs = list(pool.map(run_single_mcmc, singleChainArgs))
        
        # singleChainArgs = [i for i in tmp]
        if n_chains > 1:
            n1 = np.random.choice(range(n_chains),2,replace=False)
            [j, k] = n1
            temp_j = singleChainArgs[j][1]._temperature + 0
            temp_k = singleChainArgs[k][1]._temperature + 0
            r = (singleChainArgs[k][1]._logPost - singleChainArgs[j][1]._logPost) * temp_j + \
                (singleChainArgs[j][1]._logPost - singleChainArgs[k][1]._logPost) * temp_k
    
            # print(mc3_it, r, singleChainArgs[j][1]._logPost, singleChainArgs[k][1]._logPost, temp_j, temp_k)
            if mc3_it % 100 == 0:
                print(mc3_it, singleChainArgs[0][1]._logPost, singleChainArgs[1][1]._logPost, singleChainArgs[0][0]._w_layers[0][0][0:5])
            if r >= np.log(np.random.random()):
                singleChainArgs[j][1].reset_temperature(temp_k)
                singleChainArgs[k][1].reset_temperature(temp_j)
                print(mc3_it,"SWAPPED", singleChainArgs[j][1]._logPost, singleChainArgs[k][1]._logPost, temp_j, temp_k)

        for i in range(n_chains):
            if singleChainArgs[i][1]._temperature == 1:
                #print( singleChainArgs[i][0]._w_layers[0][0][0:10] )
                logger.log_sample(singleChainArgs[i][0],singleChainArgs[i][1])
                logger.log_weights(singleChainArgs[i][0],singleChainArgs[i][1])
            
        else:
            if mc3_it % 10 == 0:
                print(mc3_it, singleChainArgs[0][1]._logPost, singleChainArgs[0][0]._w_layers[0][0][0:5])
