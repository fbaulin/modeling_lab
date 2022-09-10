# %% Сформировать нормально распроеделенную СВ на основе равномерно распределенных СВ
"""
В программе:
1) Загрузка дополнительных библиотек
2) Инициализация параметров моделирования и распределений
3) Формирование СВ с использованем встроенной функции
4) Формировнаие СВ путём комбинирования Релеевской и арксинусной СВ
5) Задание аналитической функции распределения
6) Построение графика, содержащего:
    - аналитически заданную ПВ
    - оценённую по СВ, полученной с использованием внутренней функции
    - оценённую по СВ, полученной на основании преобразований.

Для вывода графиков в отдельное окно используйте
>>> %matplotlib qt
"""
# Загрузка дополнительных библиотек
from scipy.stats import norm, uniform, expon
import matplotlib.pyplot as plt
import numpy as np

#%%
# Моделирование СВ
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

#%% 
# Вывод ПВ результирующей нормальной СВ
# Аналитически заданная функция плотности вероятности
n_edges = 20
x_axis = np.linspace(-3*sigma_param + mu_param, 3*sigma_param + mu_param,n_edges)
w_analityc = norm.pdf(x_axis,loc=mu_param,scale=sigma_param) # аналитическое значение СВ

plt.figure()
plt.subplot(211)
plt.plot(x_axis,w_analityc,'g-',label='аналитическая')
plt.hist([y_lib,y_custom],density=True, bins=x_axis, label = ['встроенная функция','преобразование'])
plt.legend(loc='best', frameon=False)
plt.title('Нормальное распределение, m=' + str(mu_param) + ', sigma= ' + str(sigma_param))
plt.xlabel('y')
plt.ylabel('W(y)', rotation='horizontal')
plt.show()

# Вывод ПВ исходной равномерной СВ

n_edges = 40
x_axis = np.linspace(-1, 2,n_edges)
w_analityc = ((x_axis>=0) & (x_axis<=1))*1
plt.subplot(212)
plt.hlines(1,0,1, color='r', label='аналитическая ПВ')
plt.plot(x_axis,w_analityc,'g-',label='аналитическая')
plt.hist(x_1,density=True, bins=x_axis, label='исходная СВ')
plt.legend(loc='best', frameon=False)
plt.title('Гистограмма исходного процесса')
plt.xlabel('x')
plt.xlim((-1,2))
plt.ylabel('W(x)', rotation='horizontal')
plt.show()

# %% 
# Моделирование распределения заданного по варианту 
# Скопируйте текст программы в первом разделе и модифицируйте его так, чтобы получить нужное распределение

