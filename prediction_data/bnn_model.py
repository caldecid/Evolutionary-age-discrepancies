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

# set random seed
rseed = 1234
np.random.seed(rseed)


f = "/Users/dsilvestro/Software/Evolutionary-age-discrepancies/prediction_data/features.txt"
l = "/Users/dsilvestro/Software/Evolutionary-age-discrepancies/prediction_data/labels.txt"

import pandas as pd
f_tbl = pd.read_csv(f, delimiter='\t')
l_tbl = pd.read_csv(l, delimiter='\t')

#### FIT MODEL TO LOG10 AGES
# rescale features
f_tbl['Estimated_age'] = np.log10(f_tbl['Estimated_age'])
f_tbl['root_age'] = np.log10(f_tbl['root_age'])
f_tbl['number_sp'] = np.log10(f_tbl['number_sp'])
l_tbl = np.log10(l_tbl)


dat = bn.get_data(f_tbl,
                  l_tbl,
                  seed=1234,
                  testsize=0.1, # 10% test set
                  all_class_in_testset=0,
                  cv=0, # cross validation (1st batch; set to 1,2,... to run on subsequent batches)
                  header=True, # input data has a header
                  from_file=False,
                  instance_id=0,
                  randomize_order=True,
                  label_mode="regression")

raw_mse = np.mean((dat['test_labels'][:,0].flatten() - dat['test_data'][:,0].flatten())**2)

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
               n_iteration=1000000,
               sampling_f=1000,
               print_f=10000,
               n_post_samples=100,
               adapt_f=0.3,
               estimate_error=True)

# initialize output files
logger = bn.postLogger(bnn_model, filename="bnn", log_all_weights=0)

# run MCMC
bn.run_mcmc(bnn_model, mcmc, logger)


# load posterior weights
bnn_obj, mcmc_obj, logger_obj = bn.load_obj(logger._pklfile)
post_samples = logger_obj._post_weight_samples
post_weights = [post_samples[i]['weights'] for i in range(len(post_samples))]
post_alphas = [post_samples[i]['alphas'] for i in range(len(post_samples))]
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
bnn_estimates_test = np.mean(predictions_test, 0)
test_mse_bnn = np.mean((dat['test_labels'][:,0].flatten() - bnn_estimates_test)**2)

estimated_error = np.mean(np.array([post_samples[i]['error_prm'][0] for i in range(len(post_samples))]))


# run regular NN
def build_model():
    model = keras.Sequential([
        layers.Dense(8, activation='relu', input_shape=[dat['data'].shape[1]]),
        # layers.Dense(8, activation='relu'),
        layers.Dense(2, activation='relu'),
        layers.Dense(1)
    ])
    optimizer = tf.keras.optimizers.RMSprop(0.001)
    model.compile(loss='mse',
                  optimizer=optimizer,
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



# plot EXP
fig = plt.figure(figsize=(15, 5))
fig.add_subplot(131)
ax = sns.regplot(x=10**(dat['labels'][:,0].flatten()),y=10**(dat['data'][:,0].flatten()), label="Training set (Raw)")
sns.regplot(x=10**(dat['test_labels'][:,0].flatten()),y=10**(dat['test_data'][:,0].flatten()), label="Test set (MSE: %s)" % np.round(raw_mse, 3))
ax.set(xlabel='True values (Log10)', ylabel='Estimated values (Log10)')
ax.legend()
plt.axline((0, 0), (1, 1), linewidth=1, color='k')

fig.add_subplot(132)
ax = sns.regplot(x=10**(dat['labels'][:,0].flatten()),y=10**bnn_estimates_training, label="Training set (BNN)")
ax.errorbar(x=10**(dat['labels'][:,0].flatten()),y=10**bnn_estimates_training, yerr=1.96 * estimated_error, fmt='none', capsize=0, zorder=1, color='C0')
sns.regplot(x=10**(dat['test_labels'][:,0].flatten()),y=10**bnn_estimates_test, label="Test set (MSE: %s)" % np.round(test_mse_bnn, 3))
ax.errorbar(x=10**(dat['test_labels'][:,0].flatten()),y=10**bnn_estimates_test, yerr=1.96 * estimated_error, fmt='none', capsize=0, zorder=1, color='C1')
ax.set(xlabel='True values (Log10)', ylabel='Estimated values (Log10)')
ax.legend()
plt.axline((0, 0), (1, 1), linewidth=1, color='k')

fig.add_subplot(133)
ax = sns.regplot(x=10**(dat['labels'][:,0].flatten()),y=10**y_train, label="Training set (NN)")
sns.regplot(x=10**(dat['test_labels'][:,0].flatten()),y=10**y_test.flatten(), label="Test set (MSE: %s)" % np.round(test_mse, 3))
ax.set(xlabel='True values (Log10)', ylabel='Estimated values (Log10)')
ax.legend()
plt.axline((0, 0), (1, 1), linewidth=1, color='k')
fig.show()




#### FIT MODEL TO AGES (no log transformation)
# run regular NN
def build_model():
    model = keras.Sequential([
        layers.Dense(32, activation='relu', input_shape=[dat['data'].shape[1]]),
        layers.Dense(8, activation='relu'),
        layers.Dense(2, activation='relu'),
        layers.Dense(1, activation='softplus')
    ])
    optimizer = tf.keras.optimizers.RMSprop(0.001)
    model.compile(loss='mse',
                  optimizer=optimizer,
                  metrics=['mae', 'mse'])
    return model


model = build_model()
model.summary()

EPOCHS = 1000

# The patience parameter is the amount of epochs to check for improvement
early_stop = keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)

early_history = model.fit(dat['data'], 10**dat['labels'],
                          epochs=EPOCHS, validation_split=0.1, verbose=1,
                          callbacks=[early_stop])

y_train = model.predict(dat['data'], verbose=2)
y_test = model.predict(dat['test_data'], verbose=2)


test_mse = np.mean((10**(dat['test_labels'][:,0].flatten()) - y_test.flatten())**2)
test_mse_log = np.mean(((dat['test_labels'][:,0].flatten()) - np.log10(y_test).flatten())**2)

# plot results
fig = plt.figure(figsize=(10, 5))
fig.add_subplot(121)
ax = sns.regplot(x=10**dat['labels'][:,0].flatten(),y=10**dat['data'][:,0].flatten(), label="Training set (Raw)")
sns.regplot(x=10**dat['test_labels'][:,0].flatten(),y=10**dat['test_data'][:,0].flatten(), label="Test set (MSE - log: %s)" % np.round(raw_mse, 3))
ax.set(xlabel='True values', ylabel='Estimated values')
ax.legend()
plt.axline((0, 0), (1, 1), linewidth=1, color='k')
plt.xlim(0, 85)
plt.ylim(0, 85)
fig.show()


fig.add_subplot(122)
ax = sns.regplot(x=10**dat['labels'][:,0].flatten(),y=y_train, label="Training set (NN)")
sns.regplot(x=10**dat['test_labels'][:,0].flatten(),y=y_test.flatten(), label="Test set (MSE - log: %s)" % np.round(test_mse_log, 3))
ax.set(xlabel='True values', ylabel='Estimated values')
ax.legend()
plt.axline((0, 0), (1, 1), linewidth=1, color='k')
plt.xlim(0, 85)
plt.ylim(0, 85)
fig.show()
