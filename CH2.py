# Coltin Lux
# 30 November 2016 
# Mathematics/Computer Science Capstone
# 2CH Equation (Two Component Camassa - Holm Equation)
# m(sub-t) + um(sub-x) + 2mu(sub-x) = -gpp(sub-x)
# m - momentum
# u - velocity
# g - gravity
# p(row) = density

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

def findx():
    x = []
    

    for line in open("TsunamiDataRho256.txt",'r'):
        temp = line.split()
        x.append(float(temp[0]))

    return x

#==================================================================
# Find the t values that we will be testing
# Parameter is step size
#==================================================================

def findt(N):
    t = []
    index = 0
    stepsize = int(1/N)
    for i in range(0,2):
        i = float(i)
        t.append(i)
        for j in range(1,stepsize):
            index = index + N
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
    for line in open("TsunamiDataU256.txt",'r'):
        temp = line.split()
        u.append(float(temp[1]))
    return u

#==================================================================
# Initial Density
#==================================================================

def initDensity(x):
    p = []

    for line in open("TsunamiDataRho256.txt",'r'):
        temp = line.split()
        p.append(float(temp[1]))
    return p

#==================================================================
# U to M Conversion
#==================================================================

def findM(u,x,a):
    m = []
    h = (x[len(x)-1]-x[0])/len(x)

    size = len(u) - 1
    #print size

    r = []
    c = []
    temp = 1+(1/(2*(h**2)))
    temp2 = -1/(4*(h**2))
    for i in range(size):
        for j in range(size):
            if i == j:
                c.append(a*temp)
            elif ((i+2)%size) == j:
                c.append(a*temp2)
            elif ((i-2)%size) == j:
                c.append(a*temp2)
            else:
                c.append(0)
        r.append(c)
        c = []

    temp = r
    A = matrix(temp)
    # print A

    u_prime = []
    for i in range(size):
        u_prime.append(u[i+1])

    A = np.array(A)
    # print A
    u_prime = np.array(u_prime)
    # print u_prime

    C = np.dot(A,u_prime)
    #print C

    m.append(C[size-1])
    for i in C:
        m.append(i)

    return m

#==================================================================
# M to U Conversion
#==================================================================

def findU(m,x,a):
    u = []
    h = (x[len(x)-1]-x[0])/len(x)

    size = len(m) - 1
    #print size
    r = []
    c = []
    temp = 1+(1/(2*(h**2)))
    temp2 = -1/(4*(h**2))
    for i in range(size):
        for j in range(size):
            if i == j:
                c.append(a*temp)
            elif ((i+2)%size) == j:
                c.append(a*temp2)
            elif ((i-2)%size) == j:
                c.append(a*temp2)
            else:
                c.append(0)
        r.append(c)
        c = []

    temp = r
    A = matrix(temp)
    #print A

    #A inverse
    invA = A.I
    #print invA

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

def disSpace(u,m,p,g,x):
    h = (x[len(x)-1]-x[0])/len(x)
    size = len(u)-1

    m_prime = []
    p_prime = []

    for i in range(len(x)):
        if i == 0:
            firstElement = -(((u[i+1]*m[i+1])-(u[size-1]*m[size-1]))/(2*h))

            secondElement = -(m[i]*(((u[i+1]-u[size-1]))/(2*h)))

            thirdElement = ((g/2)*((((p[size-1])**2)-((p[i+1])**2))/(2*h)))

            m_prime.append(firstElement + secondElement + thirdElement)

            temp = ((((p[size-1])*(u[size-1]))-((p[i+1])*(u[i+1])))/(2*h))

            p_prime.append(temp)
        elif i == size:
            firstElement = -(((u[1]*m[1])-(u[i-1]*m[i-1]))/(2*h))

            secondElement = -(m[i]*(((u[1]-u[i-1]))/(2*h)))

            thirdElement = ((g/2)*((((p[i-1])**2)-((p[1])**2))/(2*h)))
            
            m_prime.append(firstElement + secondElement + thirdElement)

            temp = ((((p[i-1])*(u[i-1]))-((p[1])*(u[1])))/(2*h))

            p_prime.append(temp)
        else:
            firstElement = -(((u[i+1]*m[i+1])-(u[i-1]*m[i-1]))/(2*h))

            secondElement = -(m[i]*(((u[i+1]-u[i-1]))/(2*h)))

            thirdElement = ((g/2)*((((p[i-1])**2)-((p[i+1])**2))/(2*h)))
            
            m_prime.append(firstElement + secondElement + thirdElement)

            temp = ((((p[i-1])*(u[i-1]))-((p[i+1])*(u[i+1])))/(2*h))

            p_prime.append(temp)

    return m_prime,p_prime

#==================================================================
# Time Discretization
#==================================================================

def disTime(u,m,p,g,x,a):
    deltaT = .01
    a_u = []
    a_m = []
    a_p = []
    b_u = []
    b_m = []
    b_p = []
    c_u = []
    c_m = []
    c_p = []
    d_u = []
    d_m = []
    d_p = []
    e_u = []
    e_m = []
    e_p = []

    m_prime,p_prime = disSpace(u,m,p,g,x)
    u_prime = findU(m_prime,x,a)

    for i in range(len(u)):
        k = m[i] + (deltaT*m_prime[i])
        a_m.append(k)
        l = p[i] + (deltaT*p_prime[i])
        a_p.append(l)

    a_u = findU(a_m,x,a)
    
    b_m,b_p = disSpace(a_u,a_m,a_p,g,x)

    b_u = findU(b_m,x,a)

    for i in range(len(a_u)):
        k = (.75*m[i]) + (.25*a_m[i]) + (.25*deltaT*b_m[i])
        c_m.append(k)
        l = (.75*p[i]) + (.25*a_p[i]) + (.25*deltaT*b_p[i])
        c_p.append(l)

    c_u = findU(c_m,x,a)

    d_m,d_p = disSpace(c_u,c_m,c_p,g,x)

    d_u = findU(d_m,x,a)

    for i in range(len(c_u)):
        k = ((1/3)*m[i]) + ((2/3)*c_m[i]) + ((2/3)*deltaT*d_m[i])
        e_m.append(k)
        l = ((1/3)*p[i]) + ((2/3)*c_p[i]) + ((2/3)*deltaT*d_p[i])
        e_p.append(l)
    
    e_u = findU(e_m,x,a)

    return e_u,e_p

#==================================================================
# Two Component Camassa-Holm Equation (2CH)
#==================================================================

def two_CHEq(u,m,p,g,x,a):
    u_prime,p_prime = disTime(u,m,p,g,x,a)
    return u_prime,p_prime

#==================================================================
# Main Function (Initialize Variables and Call Functions)
#==================================================================

def main():
    """
    This function executes when this file is run as a script.
    """

    #==================================================================
    # Space and Time
    #==================================================================

    #space (x), we will test an array of values for x
    # 0 to 30, step size will be passed to the function
    x = findx()

    #time (t), we will test an array of values of t
    # 0 to 2, step size will be passed to the function
    t_stepsize = 0.01
    t = findt(t_stepsize)

    #==================================================================
    # Bifuration Parameter
    #==================================================================
    # smallest we can make a = .009
    a = 0.01
    
    #==================================================================
    # Find Initial Velocity
    #==================================================================
    uRecord = []
    uRecord.append(initVelocity(x))

    #==================================================================
    # Find Initial Momentum
    #================================================================== 
    mRecord = []
    mRecord.append(findM(uRecord[0],x,a))

    #==================================================================
    # Find Initial Density
    #================================================================== 
    pRecord = []
    pRecord.append(initDensity(x))

    # plt.plot(x,pRecord[0])
    # plt.show()

    #==================================================================
    # Gravity Constant
    #==================================================================

    g = 9.81

    #==================================================================
    # Call to Camassa-Holm Equation
    #==================================================================

    for i in t:
        # print (i)
        if i == 0:
            # Set to Init Values
            u = uRecord[0]
            m = mRecord[0]
            p = pRecord[0]
        else:
            u,p = two_CHEq(u,m,p,g,x,a)
            m = findM(u,x,a)
            uRecord.append(u)
            mRecord.append(m)
            pRecord.append(p)

    # print (Plotting Density)
    plt.figure(figsize=(20,15))
    plt.subplot(2,1,1)
    plt.plot(x,pRecord[2])
    plt.plot(x,pRecord[4])
    plt.plot(x,pRecord[10])
    plt.xlabel('time(s)')
    plt.ylabel('density(m/s)')
    plt.legend(('N= 2^8', 'N = 2^9','N= 2^10 '))
    plt.title('Wave Density for various N at t = 50')
 
    plt.subplot(2,1,2)
    # print (Plotting Velocity)
    plt.plot(x,uRecord[20])
    plt.plot(x,uRecord[40])
    plt.plot(x,uRecord[60])
    plt.xlabel('time(s)')
    plt.ylabel('velocity(m^3/s)')
    plt.legend(('N= 2^8', 'N = 2^9','N= 2^10 '))
    plt.title('Wave Velocity for various N at t= 50')
    plt.savefig('2CH_dp.png')
    plt.show()


# Run main() if this file is run as a script, but not if imported
if __name__ == "__main__":
    main()