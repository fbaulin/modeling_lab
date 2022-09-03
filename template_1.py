# %% Сформировать нормально распроеделенную СВ на основе равномерно распределенных СВ

from scipy.stats import norm, uniform, expon
import matplotlib.pyplot as plt
import numpy as np

# Первоначальные условия 
sigma_param = 12
mu_param=7          
size=1000
edges = 37
#------------------------------------------------------------------------------------------------------------------

#w_lib = norm.rvs(loc=mu_param,scale=sigma_param,size=size) # встроенная функция нормального распределения
w_lib = np.random.normal(loc=mu_param,scale=sigma_param,size=size) # встроенная функция нормального распределения

#x_1 = np.random.uniform(loc=0,scale=1,size=size)
x_1 = uniform.rvs(loc=0,scale=1,size=size)
x_2 = uniform.rvs(loc=0,scale=1,size=size)
y_1 = (-2 * np.log(x_1)) ** (1/2)
y_2 = np.sin(2*np.pi*x_2)
y_custom = sigma_param * y_1 * y_2 + mu_param

x = np.linspace(-3*sigma_param + mu_param, 3*sigma_param + mu_param,edges)
w_analityc = norm.pdf(x,loc=mu_param,scale=sigma_param) # аналитическое значение СВ

#-----------------------------------Вывод графиков и гистограмм-----------------------------------------------
fig, ax = plt.subplots(1, 1,  num='normilized gistogram normal function')

ax.plot(x,w_analityc,'g-',label='analytic')

ax.hist([w_lib,y_custom],density=True, bins=x, label = ['встроенная функция','custom'])
ax.legend(loc='best', frameon=False)
plt.title('нормальное распределение, m=' + str(mu_param) + ', sigma= ' + str(sigma_param))
plt.xlabel('x')
plt.ylabel('W(x)',rotation='horizontal',loc='top')
plt.show()

#-----------------------------------Вывод исходного процесса------------------------------------------------

fig, ax = plt.subplots(1, 1,  num='normilized gistogram original function')

ax.hlines(1,0,1, color='r', label='analytic')

ax.hist(x_1,density=True, bins=edges, label='original uniform')
ax.legend(loc='best', frameon=False)
plt.title('гистограмма исходного процесса')
plt.xlabel('x')
plt.ylabel('W(x)',rotation='horizontal',loc='top')
plt.show()

# %% ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-----------------------------------Моделирование заданного распределения---------------------------------------------
'''
pc_num = 4 #int(input('input PC_NUMBER: ',))

a_param = pc_num
mu_param = pc_num
sigma_param = float(pc_num / 10)
'''
lambda_param = 1
mu_param = float(1/lambda_param)
sigma_param = float(1/(lambda_param ** 0.5))

w_lib_exp = np.random.exponential(scale=mu_param, size=size) # встроенная функция экспоненциального распределения
#w_lib_exp = expon.rvs(scale=sigma_param,size=size)

X = np.linspace(0, 3 * mu_param,edges)
#w_analityc = expon.pdf(X,scale=sigma_param) # аналитическое значение СВ
w_analityc = lambda_param *  np.exp(-lambda_param * X)

x_1 = norm.rvs(loc=0,scale=1,size=size)
x_2 = norm.rvs(loc=0,scale=1,size=size)
w_custom = (x_1 ** 2 + x_2 ** 2) / (2 * lambda_param)

#--------------------------------------------------------------------------------------------------------------------
fig, ax = plt.subplots(1, 1, num='normilized gistogram exp function')

ax.plot(X,w_analityc,'g-',label='analytic')
#ax.hist(w_lib_exp,density=True, bins=edges*2, label = 'встроенная функция')
#ax.hist(w_custom,density=True, bins=edges*2, label = 'custom функция')
ax.hist([w_lib_exp,w_custom,x_1],density=True, bins=X, label = ['built-in','custom','input process'])

ax.legend(loc='best', frameon=False)
plt.title('экспоненциальная функция, лямбда= ' + str(lambda_param))
plt.xlabel('x')
plt.ylabel('W(x)', rotation='horizontal',loc='top')#va='top',ha='right')
#plt.xlim(0,3*sigma_param)
plt.show()
