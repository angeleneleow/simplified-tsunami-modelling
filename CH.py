# Camassa - Holm Equation
# m(sub-t) + um(sub-x) + 2mu(sub-x) = 0
# m - momentum
# u - velocity

# Imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numpy import matrix
from numpy import linalg

#==================================================================
# Find the x values that we will be testing
# Parameter is step size
#==================================================================

def findx(h):
    x = []
    index = 0
    stepsize = int(1/h)
    for i in range(0,30):
        i = float(i)
        x.append(i)
        for j in range(1,stepsize):
            index = index + h
            k = round((i + index), 2)
            x.append(k)
        index = 0
    x.append(30.0)
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
# Initial Velocity
#==================================================================

def initVelocity(x):
	u = []
	for i in x:
		if abs(i+5) <= 15:
			temp = (1/(np.cosh(15)))*(np.cosh(i+5))
			u.append(temp)
		if abs(i+5) > 15:
			temp = (1/(np.cosh(15)))*(np.cosh(30-(i+5)))
			u.append(temp)
	return u

#==================================================================
# Actual Solution
#==================================================================

def actualSolution(x,t):
    u = []
    for i in x:
        if abs(i+5) <= 15 + t:
            temp = (1/(np.cosh(15)))*(np.cosh(i-t+5))
            u.append(temp)
        if abs(i+5) > 15 + t:
            temp = (1/(np.cosh(15)))*(np.cosh(30-(i-t+5)))
            u.append(temp)
    return u

#==================================================================
# U to M Conversion
#==================================================================

def findM(u,x):
    m = []
    #Step Size
    h = x[len(x)-1]/len(u)

    size = len(u) - 1

    r = []
    c = []
    temp = 1+(1/(2*(h**2)))
    temp2 = -1/(4*(h**2))
    for i in range(size):
        for j in range(size):
            if i == j:
                c.append(temp)
            elif ((i+2)%size) == j:
                c.append(temp2)
            elif ((i-2)%size) == j:
                c.append(temp2)
            else:
                c.append(0)
        r.append(c)
        c = []

    temp = r
    A = matrix(temp)

    u_prime = []
    for i in range(size):
        u_prime.append(u[i+1])

    A = np.array(A)
    
    u_prime = np.array(u_prime)

    C = np.dot(A,u_prime)

    m.append(C[size-1])
    for i in C:
        m.append(i)

    return m

#==================================================================
# M to U Conversion
#==================================================================

def findU(m,x):
    u = []
    # step size
    h = x[len(x)-1]/len(m)

    size = len(m) - 1

    r = []
    c = []
    temp = 1+(1/(2*(h**2)))
    temp2 = -1/(4*(h**2))
    for i in range(size):
        for j in range(size):
            if i == j:
                c.append(temp)
            elif ((i+2)%size) == j:
                c.append(temp2)
            elif ((i-2)%size) == j:
                c.append(temp2)
            else:
                c.append(0)
        r.append(c)
        c = []

    temp = r
    A = matrix(temp)

    #A inverse
    invA = A.I

    m_prime = []
    for i in range(size):
        m_prime.append(m[i+1])

    invA = np.array(invA)
    m_prime = np.array(m_prime)

    C = np.dot(invA,m_prime)

    u.append(C[size-1])
    for i in C:
        u.append(i)

    return u

#==================================================================
# Space Discretization
#==================================================================

def disSpace(u,m,x):
    h = .1
    size = len(u)-1

    m_prime = []

    for i in range(len(x)):
        if i == 0:
            firstElement = -(((u[i+1]*m[i+1])-(u[size-1]*m[size-1]))/(2*h))

            secondElement = -(m[i]*(((u[i+1]-u[size-1]))/(2*h)))

            m_prime.append(firstElement + secondElement)
        elif i == size:
            firstElement = -(((u[1]*m[1])-(u[i-1]*m[i-1]))/(2*h))

            secondElement = -(m[i]*(((u[1]-u[i-1]))/(2*h)))
            
            m_prime.append(firstElement + secondElement)
        else:
            firstElement = -(((u[i+1]*m[i+1])-(u[i-1]*m[i-1]))/(2*h))

            secondElement = -(m[i]*(((u[i+1]-u[i-1]))/(2*h)))
            
            m_prime.append(firstElement + secondElement)

    #print m_prime
    #plt.plot(x,m_prime)
    #plt.show()
    return m_prime

#==================================================================
# Time Discretization
#==================================================================

def disTime(u,m,x):
    deltaT = .01
    a_u = []
    a_m = []
    b_u = []
    b_m = []
    c_u = []
    c_m = []
    d_u = []
    d_m = []
    e_u = []
    e_m = []

    # m and u prime are after the space discretization
    # m_prime, m
    # u_prime, u

    m_prime = disSpace(u,m,x)
    u_prime = findU(m_prime,x)

    for i in range(len(u)):
        k = m[i] + (deltaT*m_prime[i])
        a_m.append(k)

    a_u = findU(a_m,x)
    
    b_m = disSpace(a_u,a_m,x)

    b_u = findU(b_m,x)

    for i in range(len(a_u)):
        k = (.75*m[i]) + (.25*a_m[i]) + (.25*deltaT*b_m[i])
        c_m.append(k)

    c_u = findU(c_m,x)

    d_m = disSpace(c_u,c_m,x)

    d_u = findU(d_m,x)

    for i in range(len(c_u)):
        k = ((1/3)*m[i]) + ((2/3)*c_m[i]) + ((2/3)*deltaT*d_m[i])
        e_m.append(k)
    
    e_u = findU(e_m,x)
    return e_u

#==================================================================
# Camassa-Holm Equation
#==================================================================

def camassaHolmEq(u,m,x,t):
    u_prime = disTime(u,m,x)
    return u_prime

#==================================================================
# Main Function (Initialize Variables and Call Functions)
#==================================================================

def main():
    """
    This function executes when this file is run as a script.
    """

    #space (x), we will test an array of values for x
    # 0 to 30, step size will be passed to the function
    x_stepsize = .1
    x = findx(x_stepsize)

    #time (t), we will test an array of values of t
    # 0 to 2, step size will be passed to the function
    t_stepsize = .01
    t = findt(t_stepsize)
    
    #==================================================================
    # Find Initial Velocity
    #==================================================================
    uRecord = []
    uRecord.append(initVelocity(x))

    #==================================================================
    # Find Initial Momentum
    #================================================================== 
    mRecord = []
    mRecord.append(findM(uRecord[0],x))

    #==================================================================
    # Call to Camassa-Holm Equation
    #==================================================================

    for i in t:

        if i == 0:
            # Set to Init Values
            u = uRecord[0]
            m = mRecord[0]
        else:
            u = camassaHolmEq(u,m,x,t)
            uRecord.append(u)
            m = findM(u,x)
            mRecord.append(m)

    #==================================================================
    # Find Actual Solution
    #==================================================================

    uActual = []
    for i in t:
        uActual.append(actualSolution(x,i))

    #==================================================================
    # Plot
    #==================================================================

    # plt.plot(x,uRecord[0])
    # plt.plot(x,uActual[0])
    # plt.legend(('CH','actual'))
    # plt.title('Velocity plot t=0')
    # plt.show()
    # plt.plot(x,uRecord[100])
    # plt.plot(x,uActual[100])
    # plt.legend(('CH','actual'))
    # plt.title('Velocity plot at t= 100')
    # plt.show()
    # plt.plot(x,uRecord[200])
    # plt.plot(x,uActual[200])
    # plt.show()

    plt.figure(figsize = (20,15))
    plt.subplot(3,1,1)
    plt.plot(x,uRecord[0])
    plt.plot(x,uActual[0])
    plt.xlabel('time(s)')
    plt.ylabel('velocity(m^3/s)')
    plt.legend(('CH','actual'))
    plt.title('Velocity plot at t = 0')
 
    plt.subplot(3,1,2)
    plt.plot(x,uRecord[100])
    plt.plot(x,uActual[100])
    plt.xlabel('time(s)')
    plt.ylabel('velocity(m^3/s)')
    plt.legend(('CH','actual'))
    plt.title('Velocity plot at t = 100')
   
    plt.subplot(3,1,3)
    plt.plot(x,uRecord[200])
    plt.plot(x,uActual[200])
    plt.xlabel('time(s)')
    plt.ylabel('velocity(m^3/s)')
    plt.legend(('CH','actual'))
    plt.title('Velocity plot at t= 200')
    plt.savefig('CH.png')
    plt.show()


# Run main() if this file is run as a script, but not if imported
if __name__ == "__main__":
    main()