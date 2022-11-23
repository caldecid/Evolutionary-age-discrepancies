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

if __name__=="__main__":
    # set random seed
    rseed = 4321
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

    # run MCMC
    bn.run_mcmc(bnn_model, mcmc, logger)


    # load posterior weights
    bnn_obj, mcmc_obj, logger_obj = bn.load_obj(logger._pklfile)
    post_samples = logger_obj._post_weight_samples
    post_weights = [post_samples[i]['weights'] for i in range(len(post_samples))]
    post_alphas = [post_samples[i]['alphas'] for i in range(len(post_samples))]
    actFun = bnn_obj._act_fun
    output_act_fun = bnn_obj._output_act_fun

