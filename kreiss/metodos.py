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