"""
Gaussian Process Regression Implementation using Theano for symbolic gradient computation.
Author: Shen Xu
"""

# To speed Theano up, create ram disk: mount -t tmpfs -o size=512m tmpfs /mnt/randisk
# Then use flag THEANO_FLAGS=base_compiledir=/mnt/randisk python script.py
import sys, os
import theano
import theano.tensor as T
import theano.sandbox.linalg as sT
import numpy as np
import cPickle
from copy import deepcopy
import pdb

print 'Theano version: ' + theano.__version__ + ', base compile dir: ' + theano.config.base_compiledir

theano.config.mode= 'FAST_RUN'
theano.config.optimizer = 'fast_run'
theano.config.reoptimize_unpickled_function = False


def np_uniform_scalar(scale=1):
    np.random.seed(1984)
    return np.random.uniform(low=-scale,high=scale)

def shared_scalar(val=0., dtype=theano.config.floatX,name=None):
    return theano.shared(np.cast[dtype](val))


class GP_Theano(object):
    def __init__(self,
            initial_params=None):
        print 'Setting up variables ...'
        # Parameters
        if initial_params is None:
            initial_params = {'mean':None,
                              'sigma_n':0.+np_uniform_scalar(0),
                              'sigma_f':0.+np_uniform_scalar(0),
                              'l_k':0.+np.uniform_scalar(0)}
        if initial_params['mean'] == None:
            self.mean = shared_scalar(0.)
            self.meanfunc = 'zero'
        else:
            self.mean = shared_scalar(initial_params['mean'])
            self.meanfunc = 'const'
        self.sigma_n = shared_scalar(initial_params['sigma_n'])
        self.sigma_f = shared_scalar(initial_params['sigma_f'])
        self.l_k = shared_scalar(initial_params['l_k'])
        
        # Variables
        X,Y,x_test = T.dmatrices('X','Y','x_test')

        print 'Setting up model ...'
        K, Ks, Kss, y_test_mu, y_test_var, log_likelihood,L,alpha,V,fs2,sW = self.get_model(X, Y, x_test)

        print 'Compiling model ...'
        inputs = {'X': X, 'Y': Y, 'x_test': x_test}
        # solve a bug with derivative wrt inputs not in the graph
        z = 0.0*sum([T.sum(v) for v in inputs.values()])
        f = zip(['K', 'Ks', 'Kss', 'y_test_mu', 'y_test_var', 'log_likelihood',
                 'L','alpha','V','fs2','sW'],
                [K, Ks, Kss, y_test_mu, y_test_var, log_likelihood,
                 L, alpha,V,fs2,sW])
        self.f = {n: theano.function(inputs.values(), f+z, name=n, on_unused_input='ignore')
                     for n, f in f}

        if self.meanfunc == 'zero':
            wrt = {'sigma_n':self.sigma_n, 'sigma_f':self.sigma_f, 'l_k':self.l_k}
        else:
            wrt = {'mean':self.mean,'sigma_n':self.sigma_n, 'sigma_f':self.sigma_f, 'l_k':self.l_k}
        
        self.g = {vn: theano.function(inputs.values(), T.grad(log_likelihood,vv),
                                      name=vn,on_unused_input='ignore')
                                      for vn, vv in wrt.iteritems()}


    def get_model(self,X, Y, x_test):
        '''
        Gaussian Process Regression model.
        Reference: C.E. Rasmussen, "Gaussian Process for Machine Learning", MIT Press 2006

        Args:
            X: tensor matrix, training data
            Y: tensor matrix, training target
            x_test: tensor matrix, testing data
        
        Returns:
            K: prior cov matrix
            Ks: prior joint cov matrix
            Kss: prior cov matrix for testing data
            Posterior Distribution:
                alpha: alpha = inv(K)*(mu-m)
                sW: vector containing diagonal of sqrt(W)
                L: L = chol(sW*K*sW+eye(n))
            y_test_mu: predictive mean
            y_test_var: predictive variance
            fs2: predictive latent variance
        Note: the cov matrix inverse is computed through Cholesky factorization
        https://makarandtapaswi.wordpress.com/2011/07/08/cholesky-decomposition-for-matrix-inversion/
        '''
        # Compute GP prior distribution: mean and covariance matrices (eq 2.13, 2.14)
        K = self.covFunc(X,X,'K') # pior cov
        #m = T.mean(Y)*T.ones_like(Y) # pior mean
        m = self.mean*T.ones_like(Y) # pior mean

        # Compute GP joint prior distribution between training and test (eq 2.18)
        Ks = self.covFunc(X,x_test,'Ks')
        # Pay attention!! here is the self test cov matrix.
        Kss = self.covFunc(x_test,x_test,'Kss',mode='self_test')

        # Compute posterior distribution with noise: L,alpha,sW,and log_likelihood.
        sn2 = T.exp(2*self.sigma_n) # noise variance of likGauss
        L = sT.cholesky(K/sn2 + T.identity_like(K))
        sl = sn2
        alpha = T.dot(sT.matrix_inverse(L.T), 
                      T.dot(sT.matrix_inverse(L), (Y-m)) ) / sl
        sW = T.ones_like(T.sum(K,axis=1)).reshape((K.shape[0],1)) / T.sqrt(sl)
        log_likelihood = T.sum(-0.5 * (T.dot((Y-m).T, alpha)) - T.sum(T.log(T.diag(L))) - X.shape[0] / 2 * T.log(2.*np.pi*sl))
        
        
        # Compute predictive distribution using the computed posterior distribution.
        fmu = m + T.dot(Ks.T, alpha) # Prediction Mu fs|f, eq 2.25 
        V = T.dot(sT.matrix_inverse(L),T.extra_ops.repeat(sW,x_test.shape[0],axis=1)*Ks)
        fs2 = Kss - (T.sum(V*V,axis=0)).reshape((1,V.shape[1])).T # Predication Sigma, eq 2.26
        fs2 = T.maximum(fs2,0) # remove negative variance noise
        #fs2 = T.sum(fs2,axis=1) # in case x has multiple dimensions

        y_test_mu = fmu
        y_test_var = fs2 + sn2

        return K, Ks, Kss, y_test_mu, y_test_var, log_likelihood, L, alpha,V, fs2,sW

    
    def covFunc(self,x1,x2,name,method='SE',mode='cross'):
        '''
        Factorization Implementation of distance function.
        https://chrisjmccormick.wordpress.com/2014/08/22/fast-euclidean-distance-calculation-with-matlab-code/
        '''
        if method == 'SE':
            ell = T.exp(self.l_k)
            sf2 = T.exp(2.*self.sigma_f)
            if mode == 'cross':
                xx = T.sum((x1/ell)**2,axis=1).reshape((x1.shape[0],1))
                xc = T.dot((x1/ell), (x2/ell).T)
                cc = T.sum((x2/ell)**2,axis=1).reshape((1,x2.shape[0]))
                dist = xx - 2*xc + cc
            elif mode == 'self_test':
                tmp = T.sum(x1,axis=1).reshape((x1.shape[0],1))
                dist = T.zeros_like(tmp)
            else:
                raise NotImplementedError
            k = sf2 * T.exp(-dist/2)
        else:
            raise NotImplementedError
        return k

    
    def get_outputs(self, x_val, y_val, x_test_val):
        '''
        Input numpy array, output posterior distributions.
        Note: This function is independent of Theano
        '''
        inputs = {'X':x_val, 'Y':y_val, 'x_test':x_test_val}
        outputs = {n: self.f[n](*inputs.values()) for n in self.f.keys()}
        return outputs
    
    def get_prediction(self, x_val, y_val,x_test_val):
        inputs = {'X':x_val, 'Y':y_val, 'x_test':x_test_val}
        ymu = self.f['y_test_mu'](*inputs.values())
        ys2 = self.f['y_test_var'](*inputs.values())
        return ymu, ys2
    
    def get_likelihood(self,x_val, y_val):
        inputs = {'X':x_val, 'Y':y_val, 'x_test':x_val}
        likelihood = self.f['log_likelihood'](*inputs.values())
        return likelihood
 
    def get_cost_grads(self, x_val, y_val):
        '''
        get the likelihood and gradients 
        '''
        inputs = {'X':x_val, 'Y':y_val, 'x_test':x_val}
        #outputs = {n: self.f[n](*inputs.values()) for n in self.f.keys()}
        grads = {n: self.g[n](*inputs.values()) for n in self.g.keys()}
        return grads#, outputs
    

    def opt(self, train_x_val, train_y_val,params, 
            lr, momentum = 0., decay=None,
            nesterov=False, updates={},opt_method='SGD'):
        '''
        Gradient based optimizations.
        '''
        if len(updates) == 0:
            for n in params.keys():
                updates[n] = 0.
        if opt_method=='SGD':
            grads = self.get_cost_grads(train_x_val, train_y_val)
            for n in params.keys():
                g,p = grads[n], params[n]
                updates[n] = lr * g
        elif opt_method =='rmsprop':
            # RMSPROP: Tieleman, T. and Hinton, G. (2012), Lecture 6.5 - rmsprop, COURSERA:
            # Neural Networks for Machine Learning.
            if nesterov and momentum > 0.:
                # nesterov momentum, make a move according to momentum first
                # then calculate the gradients.
                for n in params.keys():
                    params[n].set_value( params[n].get_value() + momentum * updates[n])
            grads = self.get_cost_grads(train_x_val, train_y_val)
            for n in params.keys():
                g, p = grads[n], params[n]
                self.moving_mean_squared[n] = (decay * self.moving_mean_squared[n] + 
                                               (1.-decay) * g ** 2)
                updates[n] = lr * g / (np.sqrt(self.moving_mean_squared[n])+ 1e-8)
        else:
            raise NotImplementedError
        return updates


    ##############################################
    ## BEGIN TRAIN MODEL by EXTERNAL OPTIMIZERS ##
    ##############################################
    def estimate_grads(self):
        batch_size = self.batch_size
        '''
        Estimate gradient by averaging mini-batch.
        '''
        if self.meanfunc == 'zero':
            params  = {'sigma_n':self.sigma_n, 'sigma_f':self.sigma_f, 'l_k':self.l_k}
        else:
            params  = {'mean':self.mean,'sigma_n':self.sigma_n, 'sigma_f':self.sigma_f, 'l_k':self.l_k}
        N = self.x_val.shape[0]
        if batch_size is None:
            batch_size = N
        
        num_batches = N / batch_size
        if N%batch_size!=0:
            num_batches += 1
        train_index = np.arange(0,N)
        
        grads_list,est_grads = {}, {}
        for n in params.keys():
            grads_list[n] = [] 
        for i in range(num_batches):
            np.random.shuffle(train_index)
            batch_x_val = self.x_val[train_index[:batch_size],:]
            batch_y_val = self.y_val[train_index[:batch_size],:]
            grads = self.get_cost_grads(batch_x_val, batch_y_val)
            for n in params.keys():
                grads_list[n].append(grads[n])
        
        for n in params.keys():
            est_grads[n] = -1.* np.mean(grads_list[n]) # NOTE: negative grads

        return est_grads
    
    def _apply_hyp(self, hypInArray):
        '''
        Keep the order: mean, sigma_n, sigma_f, l_k
        '''
        if len(hypInArray) == 3:
            self.sigma_n.set_value(hypInArray[0])
            self.sigma_f.set_value(hypInArray[1])
            self.l_k.set_value(hypInArray[2])
        elif len(hypInArray) == 4:
            self.mean.set_value(hypInArray[0])
            self.sigma_n.set_value(hypInArray[1])
            self.sigma_f.set_value(hypInArray[2])
            self.l_k.set_value(hypInArray[3])
        else:
            raise ValueError('Number of Hyperparameters should be 3 or 4.')
    
    def _get_hypArray(self,params):
        if len(params) == 3:
            return np.array([np.sum(params['sigma_n'].get_value()),
                             np.sum(params['sigma_f'].get_value()),
                             np.sum(params['l_k'].get_value())])
        elif len(params) == 4:
            return np.array([np.sum(params['mean'].get_value()),
                             np.sum(params['sigma_n'].get_value()),
                             np.sum(params['sigma_f'].get_value()),
                             np.sum(params['l_k'].get_value())])

        else:
            raise ValueError('Number of Gradients should be 3 or 4.')

    def _convert_to_array(self, grads):
        if len(grads) == 3:
            return [grads['sigma_n'],grads['sigma_f'],grads['l_k']]
        elif len(grads) == 4:
            return [grads['mean'],grads['sigma_n'],grads['sigma_f'],grads['l_k']]
        else:
            raise ValueError('Number of Gradients should be 3 or 4.')

    def _optimizer_f(self, hypInArray):
        self._apply_hyp(hypInArray)
        ll = self.get_likelihood(self.x_val, self.y_val)
        cost = -ll # negative log-likelihood
        est_grads = self.estimate_grads()
        grads_list = self._convert_to_array(est_grads)
        return cost, np.array(grads_list)


    def train_by_optimizer(self, x_val, y_val,
            number_epoch=10, batch_size=None):
        if self.meanfunc == 'zero':
            params  = {'sigma_n':self.sigma_n, 'sigma_f':self.sigma_f, 'l_k':self.l_k}
        else:
            params  = {'mean':self.mean,'sigma_n':self.sigma_n, 'sigma_f':self.sigma_f, 'l_k':self.l_k}
        import minimize 
        if batch_size is None:
            self.batch_size = len(x_val)
        else:
            self.batch_size = batch_size
        self.x_val = x_val
        self.y_val = y_val
        
        print 'start to optimize'
        likelihood = self.get_likelihood(x_val, y_val)
        print 'BEGINE Training, Log Likelihood = %.2f'% likelihood
        
        opt_results = minimize.run(self._optimizer_f, self._get_hypArray(params),length=number_epoch,verbose=True)
        optimalHyp = deepcopy(opt_results[0])
        self._apply_hyp(optimalHyp)
        
        likelihood = self.get_likelihood(x_val, y_val)
        print 'END Training, Log Likelihood = %.2f'% likelihood
    
    ##############################################
    ## END TRAIN MODEL by EXTERNAL OPTIMIZERS   ##
    ##############################################


    def train(self, x_val, y_val,
              lr = 0.001, momentum = 0,decay = None,
              nesterov = False,batch_size=None,
              num_epoch = 10,opt_method='SGD'):
        '''
        Move hyper-parameters according to opt_method on mini-batch.
        '''
        if self.meanfunc == 'zero':
            params  = {'sigma_n':self.sigma_n, 'sigma_f':self.sigma_f, 'l_k':self.l_k}
        else:
            params  = {'mean':self.mean,'sigma_n':self.sigma_n, 'sigma_f':self.sigma_f, 'l_k':self.l_k}
        updates = {}
        # Initialize cache at the begining.
        if opt_method == 'rmsprop':
            self.moving_mean_squared={}
            for n in params.keys():
                self.moving_mean_squared[n] = 0.
        
        N = x_val.shape[0]
        if batch_size is None:
            batch_size = N
        
        num_batches = N / batch_size
        if N%batch_size!=0:
            num_batches += 1
        train_index = np.arange(0,N)
        
        #outputs = self.get_prediction(x_val, y_val, x_val) # Evaluate trained outputs
        likelihood = self.get_likelihood(x_val, y_val)
        print 'BEGINE Training, Log Likelihood = %.2f'% likelihood
        last_ll = likelihood
        for epoch in range(num_epoch):
            #if decay is not None:
            #    lr = lr * (1./(1. + decay * epoch))
            for i in range(num_batches):
                np.random.shuffle(train_index)
                batch_x_val = x_val[train_index[:batch_size],:]
                batch_y_val = y_val[train_index[:batch_size],:]

                updates = self.opt(batch_x_val, batch_y_val, params, 
                        lr=lr,momentum=momentum,decay=decay,nesterov=nesterov,updates=updates,
                        opt_method=opt_method)
                for n in params.keys():
                    p = params[n]
                    p.set_value(p.get_value() + updates[n])

            likelihood = self.get_likelihood(x_val, y_val)
            if likelihood < last_ll: # stop training
                for n in params.keys():
                    p = params[n]
                    p.set_value(p.get_value() - updates[n])
                break
            last_ll = likelihood
            print self.sigma_n.get_value(), self.sigma_f.get_value(), self.l_k.get_value()
            if epoch % 1 == 0:    
                #outputs = self.get_prediction(x_val, y_val, x_val) # Evaluate trained outputs
                print 'On Epoch %d, Log Likelihood = %.2f'%(epoch, likelihood)
        #outputs = self.get_prediction(x_val, y_val, x_val) # Evaluate trained outputs
        likelihood = self.get_likelihood(x_val, y_val)
        print 'END Training, Log Likelihood = %.2f '% likelihood






