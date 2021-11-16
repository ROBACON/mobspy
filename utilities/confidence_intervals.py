# general purpose confidence intervals, typically assumes that a NUMPY ARRAY a is passed in

import math, numpy as np
from scipy.stats import t, norm

def calc_3rd_moment(a, mean, std_dev):
	return sum([math.pow(x-mean,3) for x in a])/len(a) # 3rd sample central moment

def calc_std_dev(a, mean):
	s = math.pow(  sum([math.pow(x-mean,2) for x in a]) / (len(a)-1)  ,  1/2)
	return s

def normal_CI(a, confidence=0.95):

	mean = np.mean(a)
	s = calc_std_dev(a, mean)
	N= len(a)
	normal_conf_interval = norm.interval(confidence, loc=mean, scale=s/math.sqrt(N))

	conf_min = normal_conf_interval[0]
	conf_max = normal_conf_interval[1]
	return conf_min, conf_max


def trim_outliers(a):
	# WARNING: not really a confidence interval about sample mean
	# rather removes outlier individuals

	# rm top 5% and btm 5%, take max
	num_trim = int(len(a)*.05) #from both sides
	a = np.sort(a)
	conf_min = min(a[num_trim:-num_trim])
	conf_max = max(a[num_trim:-num_trim])
	return conf_min, conf_max


def t_CI(a, skewed=True, confidence=0.95):
	# asymmetric confidence interval based on an altered t-test
	# from eqn 2.7 of https://www.jstor.org/stable/pdf/2286597.pdf
	# thanks to https://stats.stackexchange.com/questions/16516/
	mean = np.mean(a)
	s = calc_std_dev(a, mean)
	N = len(a)
	m3 = calc_3rd_moment(a, mean, s)
	alpha = 1-confidence
	t_val = t.ppf(1-alpha/2, N-1) # alpha/2 since two-tailed

	sym_conf = t_val * s / math.pow(N,1/2)
	if skewed:
		if s == 0:
			asym_conf = 0
		else:
			asym_conf = m3 / (6*math.pow(s,2)*N)
	else:
		asym_conf = 0

	conf_min = mean - sym_conf + asym_conf
	conf_max = mean + sym_conf + asym_conf
	return conf_min, conf_max
