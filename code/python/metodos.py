import numpy as np


def choose_solver(method):
    """

    Choose solver to be used:

    Supported:

        - Euler
        - Euler improved

    """

    if method == 'euler':
        solver = euler
    elif method == 'eulerImproved':
        solver = eulerImproved
    elif method == 'RK4':
        solver = RK4

    return solver

def euler(f, k, a, b, y0):
    """

    Solve the problem

        dy/dt = f(y, t) 
        y(t=a) = y0

    in interval a <= t < b, with 
    step size k using the Euler method

    """

    time = np.arange(a, b, k)
    v_values = np.zeros(time.size)

    v = y0
    t = a
    for i, t in enumerate(time):
        v += k * f(v, t)
        v_values[i] = v

    return time, v_values

def eulerImproved(f, k, a, b, y0):
    """

    Solve the problem

        dy/dt = f(y, t) 
        y(t=a) = y0

    in interval a <= t < b, with 
    step size k using the improved 
    Euler method.

    Ref. Kreiss-Ortiz, p. 39.

    """

    time = np.arange(a, b, k)
    v_values = np.zeros(time.size)

    v = y0
    t = a
    for i, t in enumerate(time):

        v_tilde = v + k*f(v, t)
        v += k*f((v+v_tilde)/2, t+k/2) 

        v_values[i] = v

    return time, v_values   


def RK4(f, k, a, b, y0):

    time = np.arange(a, b, k)
    v_values = np.zeros(time.size)

    v = y0
    t = a
    for i, t in enumerate(time):

        q1 = f(v, t)
        q2 = f(v + (k/2)*q1, t+(k/2))
        q3 = f(v + (k/2)*q2, t+(k/2))
        q4 = f(v + k*q3, t+k)

        v += (k/6)*(q1 + 2*q2 + 2*q3 + q4) 

        v_values[i] = v

    return time, v_values   


def Qtest(f, k, a, b, y0, y_theo, method='euler'):
    
    """

    Compute precision quotient Q, defined as 

    Q(t) = ( v_1(t,k) - y(t) ) / ( v_2(t, k/2) - y(t) ).

    Ref. Kreiss-Ortiz, p. 33.

    """

    solver = choose_solver(method)
    
    time, v = solver(f, k, a, b, y0)
    time2, v2 = solver(f, k/2, a, b, y0) 
    
    num   = v - y_theo(time)
    denom = (v2 - y_theo(time2))[::2]
    
    Q = num / denom
    
    return time, Q


def QtildeTest(f, k, a, b, y0, method='euler'):
    
    """

    Compute precision quotient \tilde{Q}, defined as 

    \tilde{Q}(t) = ( v_1(t,k) - y(t) ) / ( v_2(t, k/2) - y(t) ).

    Ref. Kreiss-Ortiz, p. 33.

    """

    solver = choose_solver(method)
    
    time,  v  = solver(f, k, a, b, y0)
    time2, v2 = solver(f, k/2, a, b, y0)
    time3, v3 = solver(f, k/4, a, b, y0) 
    
    num   = v - v2[::2]
    denom = (v2 - v3[::2])[::2]
    
    Q = num / denom
    
    return time, Q

def plus(i, N):
    if i < N-1:
        return i+1
    return 0

def minus(i, N):
    if i > 0:
        return i-1
    return N-1

def Dpm_per(v, h):
    
    new_v = np.zeros_like(v)
    
    new_v[1:-1] = (v[2:] - 2*v[1:-1] + v[:-2])
    new_v[0] = v[1] - 2*v[0] + v[-1]
    new_v[-1] = v[0] - 2*v[-1] + v[-2]    

    return (1/h**2) * new_v

def D0_per(v, h):
    
    new_v = np.zeros_like(v)
    
    new_v[1:-1] = (v[2:] - v[:-2])
    new_v[0] = v[1] - v[-1]
    new_v[-1] = v[0] - v[-2]
    return (1/2*h) * new_v

def Q2_per(v, h):
    return D0_per(v, h)