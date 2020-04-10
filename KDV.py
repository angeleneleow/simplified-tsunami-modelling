# KDV Equation + Time Discretization Using Runge-Kutta
# u = u(x,t) velocity is a function of time and space
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import hypsecant

#==================================================================
# Find the x values that we will be testing
# Parameter is step size
#==================================================================

def findx(h):
    x = []
    index = 0
    stepsize = int(1/h)
    for i in range(-10,10):
        i = float(i)
        x.append(i)
        for j in range(1,stepsize):
            index = index + h
            k = round((i + index), 2)
            x.append(k)
        index = 0
    x.append(10.0)
    return x

#==================================================================
# Find the t values that we will be testing
# Parameter is step size
#==================================================================

def findt(h):
    t = []
    index = 0
    stepsize = int(1/h)
    for i in range(0,2):
        i = float(i)
        t.append(i)
        for j in range(1,stepsize):
            index = index + h
            k = round((i + index), 2)
            t.append(k)
        index = 0
    t.append(2.0)
    return t

#==================================================================
# Initial Wave Function
#==================================================================

def wavePosition(x):
    u = []
    for i in range(len(x)):
        temp = .5*((1/(np.cosh(.5*x[i])))**2)
        u.append(temp)
    return u

#==================================================================
# Actual Wave Function
#==================================================================

def actualWavePosition(x,t):
    u = []
    for i in range(len(x)):
        temp = .5*((1/(np.cosh(.5*(x[i]-t))))**2)
        u.append(temp)
    return u

#===============================================================================
# KDV Equation (to model wave movement at next time interval)
# Discretized the Ux variables
#===============================================================================

def disUX(u,x):
    # the next position
    u_prime = []

    #step size
    size = len(x)-1
    #sections
    h = 20/(len(x))

    # for each element in u, we want to calculate the derivative with respect to x
    for i in range(len(x)):
        if i == 0:
            #first calculate the discretized element 6UU(sub-x)
            firstElement = -3*((((u[i+1])**2)-((u[(size-1)])**2))/(2*h))
            #second calculate the discretized element U(sub-xxx)
            secondElement = -((u[(i+2)])-(2*(u[(i+1)]))+(2*(u[(size-1)]))-(u[(size-2)]))/((2*h)**3)
            u_prime.append(firstElement + secondElement)
        elif i == 1:
            #first calculate the discretized element 6UU(sub-x)
            firstElement = -3*((((u[i+1])**2)-((u[(i-1)])**2))/(2*h))

            #second calculate the discretized element U(sub-xxx)
            secondElement = -((u[(i+2)])-(2*(u[(i+1)]))+(2*(u[(i-1)]))-(u[(size)]))/((2*h)**3)

            u_prime.append(firstElement+secondElement)
        elif i == size:
            #first calculate the discretized element 6UU(sub-x)
            firstElement = -3*((((u[1])**2)-((u[(size-1)])**2))/(2*h))

            #second calculate the discretized element U(sub-xxx)
            secondElement = -((u[(2)])-(2*(u[(1)]))+(2*(u[(i-1)]))-(u[(i-2)]))/((2*h)**3)

            u_prime.append(firstElement+secondElement)
        elif i == (size-1):
            #first calculate the discretized element 6UU(sub-x)
            firstElement = -3*((((u[i+1])**2)-((u[(i-1)])**2))/(2*h))

            #second calculate the discretized element U(sub-xxx)
            secondElement = -((u[(1)])-(2*(u[(i+1)]))+(2*(u[(i-1)]))-(u[(i-2)]))/((2*h)**3)

            u_prime.append(firstElement+secondElement)
        else:
            #first calculate the discretized element 6UU(sub-x)
            firstElement = -3*((((u[(i+1)])**2)-((u[(i-1)])**2))/(2*h))

            #second calculate the discretized element U(sub-xxx)
            secondElement = -((u[(i+2)])-(2*(u[(i+1)]))+(2*(u[(i-1)]))-(u[(i-2)]))/((2*h)**3)

            u_prime.append(firstElement+secondElement)
    #print u_prime
    return u_prime

#===============================================================================
# KDV Equation (to model wave movement at next time interval)
# Discretized the Ut function
#===============================================================================

def disUT(u,u_ux,x):
    #change in time varible hard coded as .01
    deltaT = .005

    #list to hold all solutions that are a result of the discretized Ut function
    u_prime = u_ux

    a = []
    b = []
    c = []
    d = []
    e = []

    for i in range(len(u)):
        k = u[i] + (deltaT*u_prime[i])
        a.append(k)
    
    b = disUX(a,x)

    for i in range(len(a)):
        k = (.75*u[i]) + (.25*a[i]) + (.25*deltaT*b[i])
        c.append(k)

    d = disUX(c,x)

    for i in range(len(c)):
        k = ((1/3)*u[i]) + ((2/3)*c[i]) + ((2/3)*deltaT*d[i])
        e.append(k)

    return e


#===============================================================================
# KDV Equation (to model wave movement at next time interval
# U(sub-t) + 6UU(sub-x) + U(sub-xxx) = 0 
# However, we will not be able to easily integrate this equation
# So we will have to discretize the derivatives
#===============================================================================

def kdv(u,x):
    u_ux = disUX(u,x)
    u_ut = disUT(u,u_ux,x)
    return u_ut

#===============================================================================
# Plot Points
#===============================================================================

def plot(u,uActual,x):
    # X-Axis = x-values tested
    # Y-Axis = wave velocity or position

    plt.figure(figsize = (30,20))
    plt.subplot(3,1,1)
    plt.plot(x,u[0])
    plt.plot(x,uActual[0])
    plt.xlabel('time')
    plt.ylabel('velocity')
    plt.legend(('KDV.RK3','actual'))
    plt.title('Velocity plot at t=0')
 
    plt.subplot(3,1,2)
    plt.plot(x,u[50])
    plt.plot(x,uActual[50])
    plt.xlabel('time')
    plt.ylabel('velocity')
    plt.legend(('KDV.RK3','actual'))
    plt.title('Velocity plot at t = 100')
   
    plt.subplot(3,1,3)
    plt.plot(x,u[200])    
    plt.plot(x,uActual[200])
    plt.xlabel('time')
    plt.ylabel('velocity')
    plt.legend(('KDV.RK3','actual'))
    plt.title('Velocity plot at t= 200')
    plt.savefig('KDV_4.png')
    plt.show()

 
    return


#==================================================================
# Main Function (Initialize Variables and Call Functions)
#==================================================================

def main():
    """
    This function executes when this file is run as a script.
    """

    #space (x), we will test an array of values for x
    # -10 to 10, step size will be passed to the function
    x_stepsize = 1
    x = findx(x_stepsize)
    #print x

    #time (t), we will test an array of values of t
    # 0 to 2, step size will be passed to the function
    t_stepsize = .01
    t = findt(t_stepsize)
    #print t


    #=========================================
    # Initial Condition (time = 0)
    # Call wavePosition Funciton
    #=========================================

    #velocity (u), we will use the wave function to find the initial velocity at each x-value
    uRecord = []
    uRecord.append(wavePosition(x))

    #=========================================
    # Call to KDV Function to get position at next time interval
    #=========================================

    temp = []
    temp2 = []

    for i in range(len(t)):
        if(i == 0):
            temp.append(uRecord[0])
            temp2 = uRecord[0]
        elif(i > 0):
            temp2 = kdv(temp2,x)
            temp.append(temp2)


    u = temp

    #=========================================
    # Find the Actual Solution
    #=========================================

    uActual = []
    for i in range(len(t)):
        uActual.append(actualWavePosition(x,t[i]))

    #===============================================================================
    # Call to plot function
    #===============================================================================

    plot(u,uActual,x)




if __name__ == "__main__":
    main()