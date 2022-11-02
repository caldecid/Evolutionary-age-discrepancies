import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(suppress=True, precision=3)
# load BNN package
import sys
sys.path.insert(0, r'/Users/dsilvestro/Software/npBNN')
import np_bnn as bn
import scipy.stats
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import pandas as pd

# set random seed
rseed = 1234
np.random.seed(rseed)

# files
f =   "/Users/dsilvestro/Software/Evolutionary-age-discrepancies/test.data/ages.test.predictors.txt"
l =   "/Users/dsilvestro/Software/Evolutionary-age-discrepancies/test.data/ages.test.response.txt"
pkl = "/Users/dsilvestro/Software/Evolutionary-age-discrepancies/training.data/BNNMC3_p1_h0_l8_2_s1_binf_4594.pkl"
wd = "/Users/dsilvestro/Software/Evolutionary-age-discrepancies/test.data/"
# settings
use_estimated_rates = True
rnd_speciation_mode = False



f_tbl = pd.read_csv(f, delimiter=',')
l_tbl = pd.read_csv(l, delimiter=',')

#### FIT MODEL TO LOG10 AGES
# rescale features
f_tbl['Estimated.age'] = np.log10(f_tbl['Estimated.age'])
f_tbl['root.age'] = np.log10(f_tbl['root.age'])
f_tbl['number.sp'] = np.log10(f_tbl['number.sp'])
f_tbl['N_sisters'] = np.log10(f_tbl['N_sisters'])

# rescale labels
l_tbl = np.log10(l_tbl)


# # replace turnover and div columns with estimated ones
if use_estimated_rates:
    f_tbl['div'] = f_tbl['div.est']
    f_tbl['turnover'] = f_tbl['turnover.est']

if rnd_speciation_mode:
    f_tbl['mode'] = np.random.choice(f_tbl['mode'].to_numpy(), size=len(f_tbl['mode']), replace=False)
    f_tbl['anag'] = np.random.choice(f_tbl['anag'].to_numpy(), size=len(f_tbl['anag']), replace=False)




f_tbl_r =  f_tbl.drop(['div.est', 'turnover.est'], axis=1)



dat = bn.get_data(f_tbl_r,
                  l_tbl,
                  seed=rseed,
                  testsize=0, # 10% test set
                  all_class_in_testset=0,
                  cv=0, # cross validation (1st batch; set to 1,2,... to run on subsequent batches)
                  header=True, # input data has a header
                  from_file=False,
                  instance_id=0,
                  randomize_order=False,
                  label_mode="regression")

bnn_tmp = bn.npBNN(dat,
                     n_nodes = [8,2],
                     estimation_mode="regression",
                     actFun = bn.ActFun(fun="tanh"),
                     p_scale=1,
                    use_bias_node=3)




# load posterior weights
bnn_obj, mcmc_obj, logger_obj = bn.load_obj(pkl)


bnn_obj._test_data = bnn_tmp._data
bnn_obj._test_labels = bnn_tmp._labels


dat = {'data': bnn_obj._data,
       'labels': bnn_obj._labels,
       'test_data': bnn_obj._test_data,
       'test_labels': bnn_obj._test_labels
   }

post_samples = logger_obj._post_weight_samples
post_weights = [post_samples[i]['weights'] for i in range(len(post_samples))]
post_alphas = [post_samples[i]['alphas'] for i in range(len(post_samples))]
estimated_error = np.mean([post_samples[i]['error_prm'] for i in range(len(post_samples))])
actFun = bnn_obj._act_fun
output_act_fun = bnn_obj._output_act_fun

#### run predict
predictions_training = []
predictions_test = []
for i in range(len(post_weights)):
    actFun_i = actFun
    actFun_i.reset_prm(post_alphas[i])
    predictions_training.append(bn.RunPredict(bnn_obj._data, post_weights[i], actFun=actFun_i, output_act_fun=output_act_fun))
    predictions_test.append(bn.RunPredict(bnn_obj._test_data, post_weights[i], actFun=actFun_i, output_act_fun=output_act_fun))

predictions_training = pd.DataFrame(np.array(predictions_training)[:,:,0])
predictions_test = pd.DataFrame(np.array(predictions_test)[:,:,0])
bnn_estimates_training = np.mean(predictions_training, 0)
bnn_estimates_test = np.median(predictions_test, 0)
bnn_estimates_test_avg = np.mean(predictions_test, 0)
training_mse_bnn = np.mean((np.array(bnn_obj._labels).flatten() - bnn_estimates_training)**2)
test_labels = dat['test_labels'][:,0].flatten()
test_mse_bnn = np.mean((test_labels - bnn_estimates_test)**2)
raw_mse = np.mean((dat['test_labels'][:,0].flatten() - dat['test_data'][:,0].flatten())**2)

# calc coverage
bnn_estimates_test_M = bnn_estimates_test + 1.96 * estimated_error
bnn_estimates_test_m = bnn_estimates_test - 1.96 * estimated_error
u = np.intersect1d(np.where(test_labels[test_labels < bnn_estimates_test_M]), np.where(test_labels[test_labels > bnn_estimates_test_m]))
print("Coverage: %s" % np.round(len(u) / len(test_labels), 3) )

#### NN

#### FIT MODEL TO AGES
# run regular NN
def build_model():
    model = keras.Sequential([
        layers.Dense(32, activation='relu', input_shape=[dat['data'].shape[1]]),
        layers.Dense(8, activation='relu'),
        layers.Dense(2, activation='relu'),
        layers.Dense(1, activation='linear')
    ])
    # optimizer = tf.keras.optimizers.RMSprop(0.001)
    model.compile(loss='mse',
                  # optimizer=optimizer,
                  metrics=['mae', 'mse'])
    return model


model = build_model()
model.summary()

EPOCHS = 1000

# The patience parameter is the amount of epochs to check for improvement
early_stop = keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)

early_history = model.fit(dat['data'], dat['labels'],
                          epochs=EPOCHS, validation_split=0.1, verbose=1,
                          callbacks=[early_stop])

y_train = model.predict(dat['data'], verbose=2)
y_test = model.predict(dat['test_data'], verbose=2)
test_mse = np.mean((dat['test_labels'][:,0].flatten() - y_test.flatten())**2)
#### NN

res_list = np.array([test_labels, dat['test_data'][:,0].flatten(), bnn_estimates_test, 
            bnn_estimates_test_avg.to_numpy(), bnn_estimates_test_m.flatten(), 
            bnn_estimates_test_M.flatten(), y_test.flatten()]).T
res_tbl = pd.DataFrame(res_list)
res_tbl.columns = ["true_age", "raw_age","bnn_median","bnn_mean", "bnn_lowCI", "bnn_upCI", "nn"]
res_tbl.to_csv(os.path.join(wd, "results_bnn.csv"))

# plot results
fig = plt.figure(figsize=(15, 5))
fig.add_subplot(131)
ax = sns.regplot(x=dat['labels'][:,0].flatten(),y=dat['data'][:,0].flatten(), label="Training set (Raw)")
sns.regplot(x=dat['test_labels'][:,0].flatten(),y=dat['test_data'][:,0].flatten(), label="Test set (MSE: %s)" % np.round(raw_mse, 3))
ax.set(xlabel='True values (Log10)', ylabel='Estimated values (Log10)')
ax.legend()
plt.axline((0, 0), (1, 1), linewidth=1, color='k')
plt.xlim(-6.2, 2.5)
plt.ylim(-6.2, 2.5)

fig.add_subplot(132)
ax = sns.regplot(x=(dat['labels'][:,0].flatten()),y=bnn_estimates_training, label="Training set (BNN)")
ax.errorbar(x=(dat['labels'][:,0].flatten()),y=bnn_estimates_training, yerr=1.96 * estimated_error, fmt='none', capsize=0, zorder=1, color='C0')
sns.regplot(x=(dat['test_labels'][:,0].flatten()),y=bnn_estimates_test, label="Test set (MSE: %s)" % np.round(test_mse_bnn, 3))
ax.errorbar(x=(dat['test_labels'][:,0].flatten()),y=bnn_estimates_test, yerr=1.96 * estimated_error, fmt='none', capsize=0, zorder=1, color='C1')
ax.set(xlabel='True values (Log10)', ylabel='Estimated values (Log10)')
ax.legend()
plt.axline((0, 0), (1, 1), linewidth=1, color='k')
plt.xlim(-6.2, 2.5)
plt.ylim(-6.2, 2.5)

fig.add_subplot(133)
ax = sns.regplot(x=dat['labels'][:,0].flatten(),y=y_train, label="Training set (NN)")
sns.regplot(x=dat['test_labels'][:,0].flatten(),y=y_test.flatten(), label="Test set (MSE: %s)" % np.round(test_mse, 3))
ax.set(xlabel='True values (Log10)', ylabel='Estimated values (Log10)')
ax.legend()
plt.axline((0, 0), (1, 1), linewidth=1, color='k')
plt.xlim(-6.2, 2.5)
plt.ylim(-6.2, 2.5)
fig.show()





