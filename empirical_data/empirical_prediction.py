import os
import numpy as np
np.set_printoptions(suppress=True, precision=3)
# load BNN package
import sys
sys.path.insert(0, r'/home/torsten/Software/npBNN')
import np_bnn as bn
import pandas as pd


path = "/home/torsten/Work/Collaboration/CarlosCalderon/Evolutionary-age-discrepancies"

f = os.path.join(path, "test.data", "ages.test.predictors.txt")
f = os.path.join(path, "empirical_data", "clades_for_bnn.csv")
pkl = os.path.join(path, "training.data", "BNNMC3_p1_h0_l8_2_s1_binf_4594.pkl")
wd = os.path.join(path, "training.data")

# set random seed
rseed = 1234
np.random.seed(rseed)

f_tbl = pd.read_csv(f, delimiter = ',')

f_tbl_r = f_tbl.dropna(subset = "anag")


# load posterior weights
bnn_obj, mcmc_obj, logger_obj = bn.load_obj(pkl)


post_samples = logger_obj._post_weight_samples
post_weights = [post_samples[i]['weights'] for i in range(len(post_samples))]
post_alphas = [post_samples[i]['alphas'] for i in range(len(post_samples))]
estimated_error = np.mean([post_samples[i]['error_prm'] for i in range(len(post_samples))])
actFun = bnn_obj._act_fun
output_act_fun = bnn_obj._output_act_fun


#### run predict
predictions_mammals = np.zeros((f_tbl_r.shape[0], len(post_weights)))
for i in range(len(post_weights)):
    actFun_i = actFun
    actFun_i.reset_prm(post_alphas[i])
    predictions_mammals[:, i] = bn.RunPredict(f_tbl_r,
                                              post_weights[i],
                                              actFun = actFun_i,
                                              output_act_fun = output_act_fun).flatten()


mammals_pred_mean = np.mean(predictions_mammals, axis = 1)
mammals_pred_summary = np.stack((mammals_pred_mean,
                                 np.repeat(estimated_error, len(mammals_pred_mean)),
                                 mammals_pred_mean - 1.96 * estimated_error,
                                 mammals_pred_mean + 1.96 * estimated_error), axis = 1)
mammals_summary = pd.DataFrame(mammals_pred_summary,
                               columns = ["Estimate", "Est.error", "Q2.5", "Q97.5"])
summary_file = os.path.join(path, "empirical_data", "Mammals_PredBNN.csv")
mammals_summary.to_csv(summary_file, index = False)
