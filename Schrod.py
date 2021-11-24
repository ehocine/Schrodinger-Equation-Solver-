import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# Checking if the file already exists, so the user can overwrite it or create another one

while True:
    file_name = input('Output file name? (the file\'n format will be .dat) \n')
    file_name = file_name + '.dat'
    try:
        with open(file_name) as file:
            file_res = input('The file already exists, do you want to overwrite it? (y/n)\n')
            if file_res == 'y':
                print('The file will be overwritten...')
                break
            else:
                print('Please choose another file\'s name...')
    except FileNotFoundError:
        print('Creating file...')
        break

xmax = int(input("Max value for x (typical value: 10)? \n"))
mesh = int(input("Number of grid points (typically a few hundreds)? \n"))

dx = xmax / mesh
ddx12 = dx ** 2 / 12

# Formatting file if it already exists
data = open(file_name, 'w')
data.write('')
data.close()

data = open(file_name, 'a')  # Writing File in append mode

today_ = datetime.today().strftime('%d-%B-%Y %H:%M:%S')
data.write('File created on: ' + today_ + '\n\n')

hbar = 1
m = 1
xmu = hbar ** 2 / (2 * m)


# Potential
def V(x):
    v = 0.5 * x ** 2  # Harmonic Oscillator potential
    # v = D*(np.exp(-2*a(x-re))-2*np.exp(-a(x-re))) # Morse potential
    # v = -delta * np.exp(-delta * x) / (1 - np.exp(-delta * x)) # Hulthén potential
    # v = -1/x # Coulomb potential
    # v = -V0 / (1 + np.exp((x - R0) / a)) # Wood-Saxon potential
    # v = DD * (1 - bb * np.exp(- a * x) / (1 - np.exp(- a * x))) ** 2
    return v


# Quantum State definer

def state(n, l):
    principal_quantum_number = n + l + 1
    switcher = {
        0: "s",
        1: "p",
        2: "d",
        3: "f",
        4: "g",
        5: "h"
    }
    return str(principal_quantum_number) + switcher.get(l, "Higher State")


# Main loop
while True:

    nodes = int(input('Nodes number \n'))

    l = int(input('Orbital number \n'))

    print("State: " + state(nodes, l))

    today_ = datetime.today().strftime('%d-%B-%Y %H:%M:%S')
    data.write('Calculation started on: ' + today_ + '\n\n')
    data.write("State: " + state(nodes, l) + "\n\n")
    data.write('Nodes number: %s, Orbital number:  %s\n\n' % (nodes, l))
    data.write('Grid points number: %s\n' % mesh)

    x = list()
    pot = list()

    for i in range(0, mesh):
        x.append(i * dx)
        if x[i] == 0:
            x[i] = 1.e-4
        pot.append(V(x[i]) + l * (l + 1) * xmu / (x[i] ** 2))

    psi = np.zeros(len(x))
    f = np.zeros(len(x))

    eup = max(pot)
    elw = min(pot)

    kk = 0
    while True:

        if (eup - elw) < 1.e-10:
            break

        kk = kk + 1
        ncross = 0
        e = 0.5 * (elw + eup)
        psi[0] = x[0] ** (l + 1)
        f[0] = ddx12 * (pot[0] - e) / xmu

        icl = -1
        for i in range(1, mesh):
            f[i] = ddx12 * (pot[i] - e) / xmu

            if f[i] == 0:
                f[i] = 1.e-20
            if f[i] * f[i - 1] < 0:
                icl = i
        if icl < 0 or icl >= mesh - 2:
            eup = e
            e = 0.5 * (elw + eup)
            continue

        for i in range(len(f)):
            f[i] = 1 - f[i]

        for i in range(len(x) - 1):
            if i == 0:
                # psi[i] = (l + 1) * x[i] ** l
                psi[i + 1] = x[i + 1] ** (l + 1)
            else:
                psi[i + 1] = ((12.0 - 10.0 * f[i]) * psi[i] - f[i - 1] * psi[i - 1]) / f[i + 1]
                # psi[i + 1] = -psi[i - 1] + psi[i] * (2 - dx ** 2 * (2 * e - 2 * pot[i]))
            if psi[i] * psi[i + 1] < 0:
                ncross = ncross + 1

        print(kk, ' e = ', e)

        if ncross > nodes:
            eup = e
        else:
            elw = e

    data.write('\nIterations: %s, Bound State Energy:  %s\n\n' % (kk, e))

    data.write('Wave Function\n\n')
    data.write('x,   psi\n\n')
    for i in range(len(x)):
        data.write('%s %s\n' % (x[i], psi[i]))
    data.write('\n\n============================\n\n')
    plotting_Wave_function = input('Would you like to plot the wave function? (y/n) \n')
    if plotting_Wave_function == 'y':
        print('Plotting wave function...')
        plt.figure()
        plt.plot(x, psi)
        plt.xlabel('$r$')
        plt.xlim(0)
        plt.ylim(0)
        plt.ylabel('$\psi(r)$')
        plt.title('Wave function')
        plt.grid()
        plt.show()
    else:
        print('Not going to plot the wave function...')

    again_res = input('Would you like another calculation? (y/n)\n')

    if again_res == 'y':
        print('Launching another calculation...')
    else:
        data.close()
        print('Exiting program...')
        exit()
