import numpy as np
import random
import math
import time

def loadData(m0,sigma0,m1,sigma1,alpha0,alpha1):
    total = 1000
    da0 = np.random.normal(m0,sigma0,int(total*alpha0))
    da1 = np.random.normal(m1,sigma1,int(total*alpha1))
    dataSet = []
    dataSet.extend(da0)
    dataSet.extend(da1)
    random.shuffle(dataSet)
    return dataSet

def calcGauss(dataSetArr,mu,sigmod):
    # normal distribution pdf
    result = (1 / (math.sqrt(2*math.pi)*sigmod**2)) * np.exp(-1 * (dataSetArr-mu) * (dataSetArr-mu) / (2*sigmod**2))
    # 返回每个样本 在给定 期望，方差时的概率密度
    return result

def E_step(dataSetArr,alpha0,mu0,sigmod0,alpha1,mu1,sigmod1):
    gamma0 = alpha0*calcGauss(dataSetArr,mu0,sigmod0)
    gamma1 = alpha1*calcGauss(dataSetArr,mu1,sigmod1)
    sum_c = gamma0 + gamma1
    gamma0 = gamma0 / sum_c
    gamma1 = gamma1 / sum_c
    # 计算每个样本点 属于 每个分布的期望
    return gamma0,gamma1

def M_step(mu0, mu1, gamma0, gamma1, dataSetArr):
    mu0_new = np.dot(gamma0,dataSetArr) / np.sum(gamma0)
    mu1_new = np.dot(gamma1,dataSetArr) / np.sum(gamma1)
    
    sigmod0_new = math.sqrt(np.dot(gamma0,(dataSetArr - mu0)**2)/np.sum(gamma0))
    sigmod1_new = math.sqrt(np.dot(gamma1,(dataSetArr - mu1)**2)/np.sum(gamma1))
    alpha0_new = np.sum(gamma0) / len(gamma0)
    alpha1_new = np.sum(gamma1) / len(gamma1)
    return mu0_new,mu1_new,sigmod0_new,sigmod1_new,alpha0_new,alpha1_new

def EM_Train(dataSetList,iter=1000):
    dataSetArr = np.array(dataSetList)
    # 迭代结果与初始值设置密切相关
    alpha0 = 0.5
    mu0 = -1
    sigmod0 = 1
    alpha1 = 0.5
    mu1 = 1
    sigmod1 = 1
    step = 0
    while(step < iter):
        gamma0, gamma1 = E_step(dataSetArr, alpha0, mu0, sigmod0, alpha1, mu1, sigmod1)
        mu0, mu1, sigmod0, sigmod1, alpha0, alpha1 = M_step(mu0, mu1, gamma0, gamma1, dataSetArr)
        step += 1
    return alpha0, mu0, sigmod0, alpha1, mu1, sigmod1

start = time.time()
alpha0 = 0.3  
mu0 = -2  
sigmod0 = 0.5  
alpha1 = 0.7  
mu1 = 3  
sigmod1 = 1  
dataSetList = loadData(mu0, sigmod0, mu1, sigmod1, alpha0, alpha1)
print('alpha0:%.1f, alpha1:%.1f, mu0:%.1f, mu1:%.1f, sigmod0:%.1f, sigmod1:%.1f' % (alpha0, alpha1, mu0, mu1, sigmod0, sigmod1 ))
alpha0, mu0, sigmod0, alpha1, mu1, sigmod1 = EM_Train(dataSetList)
print("---------------")
print('alpha0:%.1f, alpha1:%.1f, mu0:%.1f, mu1:%.1f, sigmod0:%.1f, sigmod1:%.1f' % (alpha0, alpha1, mu0, mu1, sigmod0, sigmod1))
print('time span:', time.time() - start)

# example 2
from scipy import stats
import numpy as np
# 硬币投掷结果观测序列,求解 两个硬币投出正面的概率
observations = np.array([[1, 0, 0, 0, 1, 1, 0, 1, 0, 1],
                         [1, 1, 1, 1, 0, 1, 1, 1, 1, 1],
                         [1, 0, 1, 1, 1, 1, 1, 0, 1, 1],
                         [1, 0, 1, 0, 0, 0, 1, 1, 0, 0],
                         [0, 1, 1, 1, 0, 1, 1, 1, 0, 1]])
def em_single(priors, observations):

    counts = {'A': {'H': 0, 'T': 0}, 'B': {'H': 0, 'T': 0}}
    theta_A = priors[0]
    theta_B = priors[1]
    # E step
    for observation in observations:
        len_observation = len(observation)
        num_heads = observation.sum()
        num_tails = len_observation - num_heads
        contribution_A = stats.binom.pmf(num_heads, len_observation, theta_A)
        contribution_B = stats.binom.pmf(num_heads, len_observation, theta_B)   # 两个二项分布
        weight_A = contribution_A / (contribution_A + contribution_B)
        weight_B = contribution_B / (contribution_A + contribution_B)
        # 更新在当前参数下A、B硬币产生的正反面次数
        counts['A']['H'] += weight_A * num_heads
        counts['A']['T'] += weight_A * num_tails
        counts['B']['H'] += weight_B * num_heads
        counts['B']['T'] += weight_B * num_tails
    # M step
    new_theta_A = counts['A']['H'] / (counts['A']['H'] + counts['A']['T'])
    new_theta_B = counts['B']['H'] / (counts['B']['H'] + counts['B']['T'])
    return [new_theta_A, new_theta_B]

def em(observations, prior, tol=1e-6, iterations=10000):
    iteration = 0
    while iteration < iterations:
        new_prior = em_single(prior, observations)
        delta_change = np.abs(prior[0] - new_prior[0])
        if delta_change < tol:
            break
        else:
            prior = new_prior
            iteration += 1
    return [new_prior, iteration]
 
if __name__ == "__main__":
    result=em(observations, [0.6, 0.4])
    print(result)
