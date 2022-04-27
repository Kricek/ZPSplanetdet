import numpy as np

A = 0.61
q_br = 1.7 * 10**(-4)
n = -0.93
p = 0.6
m = 0.49

def calculate_f(s, q, params):
    """
    this function calculates f(s, q, A, q_br, n, p, m) = A*s**m*[(q/q_br)**n*theta(q-q_br)+(q/q_br)**p*theta(q_br-q)]
    where theta is the Heaviside function
    """
    if(q > params[1]):
        f = params[0] * ((q/params[1]) ** (params[2])) * (s ** params[4])
    else:
        f = params[0] * ((q/params[1]) ** (params[3])) * (s ** params[4])
    return f

params = (A, q_br, n, p, m)
(q, s) = np.loadtxt("sq_dgen.dat", usecols = (1, 2), unpack = True)

f = np.zeros((len(q), 1))
ratio = np.zeros((int(len(q) / 2), 1))

k = 0
	
for i in range(0, len(q)):	
    f[i] = calculate_f(s[i], q[i], params)

for j in range(0, int(len(q) / 2)):
    ratio[j] = f[k] / f[k+1]
    print(*ratio[j])
    k += 2
	
