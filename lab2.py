import numpy as np
from scipy.integrate import odeint
import matplotlib.pylab as plt

# graphs for differential equations of question1 and question2 of lab2
# x_q1, x_q2 are the solutions in form x(t)=C1cos(nt)+C2*sin(nt)+f(t) where C1 =c1+c2 and C2=i(c1-c2)
# q1: x''(t)+(k/m)*x(t)=g+(1/m)*F*cos(w*t)
# q2: x''(t)+(k/m)*x(t)=g+(1/m)*(k/m)^(1/2)

def x_q1(t,x0,xpi2):
    f_pi2 = (g/(n**2))+(F*np.sin((w*np.pi)/2)*np.sin((n/2)*np.pi))/(2*n*m*w)

    #constant variables of function, for C1 = C2* and C2 = a+bi, where a,b real number 
    a = (1/2)*(x0-(g/(n**2))) # x0: initial condition for t=0
    b = (1/2)*(xpi2-f_pi2) # xpi2: initial condition for t=pi/2

    f = (g/(n**2))+ (F*np.sin(w*t)*np.sin(n*t))/(2*n*m*w)
    return 2*a*np.cos(n*t)+2*b*np.sin(n*t)+f


def x_q2(t,x0,xpi2):
    f = (m*g + n)/k

    #constant variables of function, for C1 = C2* and C2 = a+bi, where a,b real number 
    a = (1/2)*(x0 - f) # x0: initial condition for t=0
    b = (1/(2*np.sin((n/2)*np.pi)))*(xpi2 - f) # xpi2: initial condition for t=pi/2

    return 2*a*np.cos(n*t)+2*b*np.sin(n*t)+f


def makeSubplot(l, x0_1, xpi2_1, x0_2, xpi2_2, gr_title): # l: is the height of the tuple in subplot which define the starting point of the top left corner
    sub_plt = plt.subplot2grid((11,1), (l,0), rowspan = 3, colspan = 1)
    
    if (l == 8): # 3rd subplot
        sub_plt.set_xlabel('t',fontsize=14)

    soln_q1 = x_q1(time,x0_1,xpi2_1)
    soln_q2 = x_q2(time,x0_2,xpi2_2)

    sub_plt.set_title(gr_title)
    sub_plt.plot(time, soln_q1, label = 'Q1', color = 'b')
    sub_plt.plot(time, soln_q2, label = 'Q2', color = 'g')
    sub_plt.legend()


def makeGraph(x0_1, x0_2, graph_title):
    # initial condition for t=pi/2
    # - - - eq_of_question1 - - - #
    f_pi2 = (g/(n**2))+(F*np.sin((w*np.pi)/2)*np.sin((n/2)*np.pi))/(2*n*m*w)
    xpi2_pos1 = 6*f_pi2             # C2 > 0, this happens for xpi2>g/(n**2)
    xpi2_01   = f_pi2               # C2 = 0, this happens for xpi2=g/(n**2)
    xpi2_neg1 = (1/3)*f_pi2         # C2 < 0, this happens for xxpi2<g/(n**2)
    # - - - eq_of_question2 - - - #
    xpi2_pos2 = 10*((m*g)+n)/k      # C2 > 0, this happens for xpi2>((m*g)+n)/k
    xpi2_02   = (m*g + n)/k         # C2 = 0, this happens for xpi2=((m*g)+n)/k
    xpi2_neg2 =(1/3)*((m*g)+n)/k    # C2 < 0, this happens for xpi2<((m*g)+n)/k

    # - C2>0 (1st subplot)
    makeSubplot(0, x0_1, xpi2_pos1, x0_2, xpi2_pos2, graph_title + 'C2>0')
    # - C2=0 (2nd subplot)
    makeSubplot(4, x0_1, xpi2_01, x0_2, xpi2_02, graph_title + 'C2=0')
    # - C2<0 (3rd subplot)
    makeSubplot(8, x0_1, xpi2_neg1, x0_2, xpi2_neg2, graph_title + 'C2<0')

    plt.show() 


# - - - Main program - - - #
#parameters of problem
m = 10
k = 10**2
n = (k/m)**(1/2)
w = 5
g = 9.81
F = 10


# initial condition (where C1=c1+c2, C2=c1-c2 based on-hand solved differential equations)
# - - - equation of question1 - - - #
x0_pos1 = 3*g/(n**2)            # C1 > 0, this happens for x0>g/(n**2)
x0_01 = g/(n**2)                # C1 = 0, this happens for x0=g/(n**2)
x0_neg1 =(1/10)*g/(n**2)        # C1 < 0, this happens for x0<g/(n**2)

# - - - equation of question2 - - - #
x0_pos2 = 3*((m*g)+n)/k         # C1 > 0, this happens for x0>((m*g)+n)/k
x0_02 = (m*g + n)/k             # C1 = 0, this happens for x0=((m*g)+n)/k
x0_neg2 =(1/30)*((m*g)+n)/k     # C1 < 0, this happens for x0<((m*g)+n)/k


# perform numerical integration
# time instances for calculating solution
time_i = 0 
time_f = 3.5
time = np.arange(time_i,time_f,0.2)

# plot the numerical solution
# - for C1>0 (1st graph)
makeGraph(x0_pos1, x0_pos2, 'C1>0 and ')
# - for C1=0 (2nd graph)
makeGraph(x0_01, x0_02, 'C1=0 and ')
# - for C1<0 (3rd graph)
makeGraph(x0_neg1, x0_neg2, 'C1<0 and ')