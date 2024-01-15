import matplotlib.pyplot as plt

t, ekin, poteng, toteng = [], [], [], []

with open('thermo.log', 'r+') as f1:
    lines = f1.readlines()
    for i,line in enumerate(lines):
        if i > 0:
            row = line.split()
            t.append(float(row[1]))
            ekin.append(float(row[4]))
            poteng.append(float(row[3])) 
            toteng.append(float(row[5]))

fig1 = plt.figure()
plt.title('System energies as a function of time')
plt.plot(t, toteng, label='Total Energy')
plt.plot(t, ekin, label='Kinetic Energy')
plt.plot(t, poteng, label='Potential Energy')
plt.xlabel('Time (fs)')
plt.ylabel('Energy (eV)')
plt.legend()
plt.savefig('Energy.png', dpi=100)