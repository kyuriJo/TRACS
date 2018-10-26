
import numpy as np
import matplotlib
matplotlib.use("Agg")
import gapGP
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern,RBF,WhiteKernel
from sklearn.metrics import mean_squared_error

def multipleReg(X, X_pred, Ys, pjname, subname, c, i, minmax, bound=(0.1, 0.1)) :
  #print "Start regression of ", pjname, subname, c, i
  lbound = (1.0, np.max(X))
  X_ = X.reshape((-1,1))
  X_pred_ = np.array(X_pred).reshape((-1,1))
  kern = 1.0 * RBF(length_scale=lbound[0], length_scale_bounds=lbound)+WhiteKernel(noise_level=bound[0], noise_level_bounds=bound)
  Ym = np.mean(Ys, axis=0)
  Yv = np.var(Ys, axis=0)
  m = GaussianProcessRegressor(kernel=kern, n_restarts_optimizer=10)
  m.fit(X_, Ym)
  p = m.kernel_.get_params()
  params = [p['k1__k1__constant_value'], p['k1__k2__length_scale'], m.log_marginal_likelihood(m.kernel_.theta)]
  like = m.log_marginal_likelihood(m.kernel_.theta)
  pred_m, pred_v = m.predict(X_pred_, return_std=True)
  pred_m_X, pred_v_X = m.predict(X_, return_std=True)
  if (subname!='TA'):
    fig = plt.figure()
    mse = mean_squared_error(Ys, np.repeat(pred_m_X.reshape((1, len(X))), Ys.shape[0], axis=0))
    dist = np.mean([gapGP.calDist_l(Ys[ii,:], X, X_pred, pred_m, pred_v) for ii in range(Ys.shape[0])])
    for y in Ys :
      plt.plot(X_, y, linestyle='solid', color='0.5', zorder=1)
    plt.plot(X_pred, pred_m, 'r-', zorder=9)
    plt.plot(X_, Ym, 'r.', zorder=10)
    plt.fill(np.concatenate([X_pred, X_pred[::-1]]),
       np.concatenate([pred_m - 1.9600 * pred_v, (pred_m + 1.9600 * pred_v)[::-1]]),
       alpha=.3, fc='r', ec='None', zorder=8)
    plt.xlabel('$x$')
    plt.ylabel('$f(x)$')
    plt.ylim(minmax[0], minmax[1])
    plt.title("Posterior (kernel: %s)\n Log-Likelihood: %.3f , avgMSE: %.3f, avgLike: %.3f"
            % (m.kernel_, like, mse, dist),
            fontsize=12)
    plt.savefig(pjname+'/'+subname+'P'+str(c)+'cluster'+str(i)+'.png')
    plt.close()
  return pred_m_X, pred_v_X, pred_m, pred_v, params

def regression(X, X_pred, Y, bound):
  kern = 1.0 * Matern(length_scale=1.0, length_scale_bounds=(1e-1, 720), nu=1.5)
  kern = 1.0 * RBF(length_scale=100.0, length_scale_bounds=(5,100))+WhiteKernel(noise_level=bound[0], noise_level_bounds=bound)
  X_ = np.array(X).reshape((-1,1))
  X_pred_ = np.array(X_pred).reshape((-1,1))
  m = GaussianProcessRegressor(kernel=kern, n_restarts_optimizer=10)
  m.fit(X_, Y)
  pred_m, pred_v = m.predict(X_pred_, return_std=True)
  return pred_m, pred_v
