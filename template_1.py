# %% Сформировать нормально распроеделенную СВ на основе равномерно распределенных СВ
"""
В программе:
1) загрузка дополнительных библиотек
2) Инициализация параметров моделирования и распределений
3) Формирование СВ с использованем встроенной функции
4) Формировнаие СВ путём комбинирования Релеевской и арксинусной СВ
5) Задание аналитической функции распределения
6) Построение графика, содержащего:
    - аналитически заданную ПВ
    - оценённую по СВ, полученной с использованием внутренней функции
    - оценённую по СВ, полученной на основании преобразований.
"""
from scipy.stats import norm, uniform, expon
import matplotlib.pyplot as plt
import numpy as np

# Параметры моделирования и распределения
sigma_param = 12
mu_param=7          
n_values=1000

# Формирование нормально распределенных СВ
y_lib = np.random.normal(loc=mu_param,scale=sigma_param,size=n_values) # встроенная функция нормального распределения
# Формиронование нормально распределенной СВ (см. раздел 1.2.7)
x_1 = np.random.uniform(0,1,size=n_values)      # исходные равномерные СВ  
x_2 = np.random.uniform(0,1,size=n_values)      # исходные равномерные СВ
y_1 = (-2 * np.log(x_1)) ** (1/2)       # результат н.лин. преобразования - Релеевская СВ
y_2 = np.sin(2*np.pi*x_2)               # результат н.лин. преобразования - Арксинусная СВ
y_custom = sigma_param * y_1 * y_2 + mu_param   # комбинирование для получения нормальной СВ

# Аналитически заданная функция плотности вероятности
n_edges = 20
x = np.linspace(-3*sigma_param + mu_param, 3*sigma_param + mu_param,n_edges)
w_analityc = norm.pdf(x,loc=mu_param,scale=sigma_param) # аналитическое значение СВ

# Вывод ПВ результирующей нормальной СВ
fig, ax = plt.subplots(1, 1,  num='normilized histogram normal function')
ax.plot(x,w_analityc,'g-',label='аналитическая')
ax.hist([y_lib,y_custom],density=True, bins=x, label = ['встроенная функция','преобразование'])
ax.legend(loc='best', frameon=False)
plt.title('нормальное распределение, m=' + str(mu_param) + ', sigma= ' + str(sigma_param))
plt.xlabel('y')
plt.ylabel('W(y)',rotation='horizontal',loc='top')
plt.show()

# Вывод ПВ исходной равномерной СВ
n_edges = 40
x = np.linspace(-1, 2,n_edges)
fig, ax = plt.subplots(1, 1,  num='normilized histogram original function')
ax.hlines(1,0,1, color='r', label='analytic')
ax.hist(x_1,density=True, bins=x, label='original uniform')
ax.legend(loc='best', frameon=False)
plt.title('гистограмма исходного процесса')
plt.xlabel('x')
plt.xlim((-1,2))
plt.ylabel('W(x)',rotation='horizontal',loc='top')
plt.show()

# %% 
# Моделирование распределения заданного по варианту 
# Скопируйте текст программы в первом разделе и модифицируйте его так, чтобы получить нужное распределение

