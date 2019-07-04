from model import Model

params0 = params
model0 = Model(params0)
model0.build()

t = [0.5, 2.5, 7.5,9.9]
model0.plot(t,'E')
model0.plot(t,'phi')

#Compute surface electric field in terms of the bulk concentration of electrolytes
numPoints = 50
Cb = np.linspace(0, 1.4, numPoints)
surfE_C = []
aux1 = []
aux2 = []
aux3 = []
t = [0.1, 0.5, 0.7]
params1 = params
for i in range(len(Cb)):
    C = i
    params1['bulkConcentration'] = C
    model3 = Model(params1)
    model3.build()
 
    n = int(t[0]/model3.dtau)
    aux1.append(model3.E[n,0])
    n = int(t[1]/model3.dtau)
    aux2.append(model3.E[n,0])
    n = int(t[2]/model3.dtau)
    aux3.append(model3.E[n,0])
    
surfE_C.append([aux1, aux2, aux3])

surfE_C =np.array(surfE_C)/1e6
surfE_C = np.array(surfE_C)

#Compute surface electric field in terms of the electrode's voltage
surfE_v = []
aux1 = []
aux2 = []
aux3 = []
v = np.linspace(-1.5, 0, numPoints)
params1['bulkConcentration'] = 1e-3
for i in range(0,len(v)):
    V0 = v[i]
    params1["V0"] = V0
    model4 = Model(params1)
    model4.build()
    n = int(t[0]/model4.dtau)
    aux1.append(model4.E[n,0])
    n = int(t[1]/model4.dtau)
    aux2.append(model4.E[n,0])
    n = int(t[2]/model4.dtau)
    aux3.append(model4.E[n,0])
    
surfE_v.append([aux1, aux2, aux3])
surfE_v = np.array(surfE_v)

with open('surfE_v.txt') as v_file:
    for value in surfE_v:
        v_file.write(value)
with open('surfE_v.txt') as c_file:
    for value in surfE_C:
        c_file.write(value)
        
#Plot results
fig, ax1 = plt.subplots()

mw = 4
fs = 26
skip = 4
t = [0.1, 0.5, 0.7]
plt.style.use('thesis')
plt.title('Surface Electric Field As A Function Of Concentration')
plt.legend(loc='upper left')
for i in range(len(t)):

    plt.plot(Cb, surfE_C[0][i], '--', label=r'$E_s$,  $t =%.2f ns$' % Decimal(t[i] * 1e9/ ( model4.D1 * model4.kappa ** 2 ) ))
    
plt.xlabel(r'Molar Concentration', fontsize=fs+5)
plt.ylabel(r'Surface Electric Field (V/m) $\times 10^6$')
textstr = '\n'.join((
                r'$k_f=%.2f m/s$' % Decimal(model3.kf),
                r'$\kappa=%.2f nm^{-1}$' % Decimal(model3.kappa/1e9),
                r'$V_0=%.2f V$' % Decimal(model3.V0),))

props = dict(boxstyle='round', facecolor='white', alpha=0.5)
plt.text(0.6, 0.95, textstr, transform=ax1.transAxes, fontsize=fs, verticalalignment='top', bbox=props)
plt.legend(loc='upper right', fontsize = fs)
plt.xticks(size = fs)
plt.yticks(size = fs)
plt.savefig('../../img/surfaceEfield_Cb.eps', dpi=1000, fontweight='bold')
#plt.show()

plt.figure()
plt.style.use('thesis')
skip = 4
plt.title('Surface Electric Field As A Function Of V', fontweight='bold')
plt.legend(loc='upper left')
for i in range(len(t)):
    n = int(t[i]/model3.dtau)
    plt.plot(v, surfE_v[0][i], label=r'$E_s$,  $t =%.2f ns$' % Decimal(t[i] * 1e9/ ( model4.D1 * model4.kappa ** 2 ) ))
plt.xlabel(r'Surface Voltage (V)')
plt.ylabel(r'Surface Electric Field (V/m) $\times 10^8$')
textstr = '\n'.join((
                r'$k_f=%.2f$' % Decimal(model4.kf),
                r'$\kappa=%.2f m^{-1}$' % Decimal(model4.kappa),
                r'$C_b=%.2f M$' % Decimal(1e3 * model4.Cb),))

props = dict(boxstyle='round', facecolor='white', alpha=0.5)
plt.text(0.68, 0.3, textstr, transform=ax1.transAxes, fontsize=30, verticalalignment='top', bbox=props)

plt.legend(loc='lower right')
plt.savefig('../../img/surfaceEfield_v.eps', dpi=1000, fontweight='bold')
plt.show()


Cbulk = 1e-3*np.array([0.01, 0.3, 0.7, 1.2])
deltaE = []
timewindows = []
#new time span is an array which contains the wanted time span in ns 
#(in this case 50 ns) and then transforms it to adimensional units via
# the factor D * kappa. This must be calculated first and the come back
#to this part of the code to reinsert the kappas.
new_time_span = 50 * 1e-9 * np.array([1/(2.0382813587050528 * 1e-10), 1/(2.9894793261007435* 1e-10),1/(4.484218989151115 * 1e-10)])
#new_time_span = 50 * 1e-9 * np.array([1/(4.484218989151115 * 1e-10), 2/(3 * 2.9894793261007435* 1e-10), 1/ (2.2 * 2.0382813587050528 * 1e-10)])
params2 = params
for C in Cbulk:
    params2['bulkConcentration'] = C
    kappa = np.sqrt(( params2["z"] * params2["Fa"] ) ** 2 * params2["bulkConcentration"] / ( params2["epsilon"] * params2["T"] * params2["R"] ))
    params2["timespan"] = ( params2["D1"] * params2["diffusionCoefficientScale"] * kappa ** 2 ) * 50 * 1e-9

    model5 = Model(params2)
    model5.build()
    timewindows.append( params2["timespan"] / ( model5.D1 * params2["diffusionCoefficientScale"] * model5.kappa ** 2 ) )
    deltaE.append(model5.E[:] - model5.E[-1])# * np.zeros(len(model5.E[1:,0])))


deltaE = np.array(deltaE)#/1e6


t = []
for i in range(len(Cbulk)):
     t.append(np.linspace(0, timewindows[i], len(deltaE[0])))


        
plt.figure(figsize=(20,16))
plt.grid(True, color= '#F2F2F2')
mw = 4
fs = 26
skip = 4

for i in range(len(Cbulk)):
    plt.plot(t[i][1:], deltaE[i][1:,0], '-', label=r'$E_s$,  $C_b ='+str(1e3*Cbulk[i])+'M$')
    
plt.xlabel(r'Time (ns)')
plt.ylabel(r'Electric Field Fluctuations (V/m) $\times 10^6$')
plt.title(r'Electric Field Fluctuations at the electrode $\delta E$')
textstr = '\n'.join((
                r'$k_f=%.2f$' % Decimal(model4.kf),
                r'$\kappa=%.2f m^{-1}$' % Decimal(model4.kappa),
                r'$V_0=%.2f M$' % Decimal(model4.V0),))

props = dict(boxstyle='round', facecolor='white', alpha=0.5)
plt.text(0.68, 0.7, textstr, transform=ax1.transAxes, fontsize=30, verticalalignment='top', bbox=props)

plt.legend(loc='upper right')
plt.savefig('../../img/surfaceDeltaE.eps', dpi=1000, fontweight='bold')
plt.show()


#Fit Curve
from scipy.optimize import curve_fit 

def polinomial(x, a, c):
    b=1
    return a * x ** b + c

def exponencial(x, a, b):
    return a * np.exp(-b * x)

fitParams = []
pcovs = []
for i in range(len(t)):
    fit_params, pcov = curve_fit(polinomial, t[i][1:], deltaE[i][1:,0])
    fitParams.append(fit_params)
    pcovs.append(pcov)


plt.figure(figsize=(20,16))

plt.title('Electric Field Fluctiation Polinomial Fit', fontweight='bold')
plt.legend(loc='upper left')



for i in range(len(Cbulk)):
    a = fitParams[i][0]
    b = fitParams[i][1]
    plt.plot(t[i][1:], deltaE[i][1:,0], '-', label=r'$E_s$,  $C_b ='+str(1e3*Cbulk[i])+'M$')
    plt.plot(t[i][1:], polinomial(t[i][1:],a,b) , '--')


plt.xlabel(r'Time (ns)')
plt.ylabel(r'Electric Field Fluctuations (V/m)')

textstr = '\n'.join((
                r'$k_f=%.2f$' % Decimal(model4.kf),
                r'$\kappa=%.2f m^{-1}$' % Decimal(model4.kappa),
                r'$V_0=%.2f M$' % Decimal(model4.V0),
                r'$C_b=%.2f M$' % Decimal(model4.Cb),),)

#props = dict(boxstyle='round', facecolor='white', alpha=0.5)
#plt.text(0.65, 0.1, textstr, transform=ax1.transAxes, fontsize=fs-4, verticalalignment='top', bbox=props)

plt.legend(loc='upper right')
plt.savefig('../../img/surfaceDeltaE_fitting.eps', dpi=1000, fontweight='bold')
plt.show()
print(fitParams)