class Model():
    def __init__(self, params):
        #Model Parameters
        self.Cb = params['bulkConcentration']
        self.D = params['diffusionCoefficientScale']
        self.d = params['laminarFlowRegion']
        self.kf = params['reactionRate']
        self.z = params['z']
        self.F = params['Fa']
        self.R = params['R']
        self.T = params['T']
        self.V0 = params['V0']
        self.epsilon = params['epsilon']
        self.kappa =  np.sqrt(( ( self.z * self.F  ) ** 2 * self.Cb ) / ( self.epsilon * self.R * self.T ) )
        self.Psi0 = self.z * self.F * params['V0'] / ( self.R * self.T )
        self.D1 = self.D * params["D1"]
        self.D2 = self.D * params["D2"]
        self.N = 30000
        self.M = 100
        self.length = params["length"]
        self.xi = np.linspace(0,params["length"], self.M)
        self.tau = np.linspace(0,params["timespan"], self.N) #shape is N+1

        #Grid Parameters
        self.dtau = params["timespan"]/(self.N)  # N Partitions
        self.dxi = self.length/(self.M) # M Partitions 
        self.a1 = self.dtau / self.dxi ** 2 
        self.a2 = self.dtau / self.dxi ** 2 * self.D2/self.D1

        #Plotting parameters

    def build(self):
        M = self.M
        N = self.N
        a1 = self.a1
        a2 = self.a2
        Psi0 = self.Psi0
        kappa = self.kappa
        kf = self.kf
        dxi = self.dxi
        D1 = self.D1 
        D2 = self.D2
        
        # Define the coefficient matrix
        g1 = 1 / ( 1 + kf * dxi / ( D1 * kappa ) +  Psi0)
        di1 = ( 1 - 2 * a1 ) * np.ones(M-2)
        di1[0] = ( 1 - 2 * a1 + a1 * g1 )
        A1 = diags(np.array([ a1 * np.ones(M-3), di1, a1 * np.ones(M-3)]), [-1, 0, 1], shape=(M-2, M-2)).toarray()

        g2 = 1 / ( 1 - Psi0 )
        di2 = ( 1 - 2 * a2 ) * np.ones(M-2)
        di2[0] = ( 1 - 2 * a2 + a2 * g2 )
        A2 = diags(np.array([ a2 * np.ones(M-3), di2, a2 * np.ones(M-3)]), [-1, 0, 1], shape=(M-2, M-2)).toarray()

        B1 = np.zeros([M-2, M-2])
        B2 = np.zeros([M-2, M-2])

        D0 = diags(np.array([ np.ones(M-3), -2 * np.ones(M-2), np.ones(M-3)]), [-1, 0, 1], shape=(M-2, M-2)).toarray()
        Dinv = np.asarray(np.linalg.inv(D0))

        b1 = np.zeros(M-2)
        b1[-1] = a1 
        b2 = np.zeros(M-2)
        b2[-1] = a2

        bPsi = np.zeros(M-2)
        bPsi[0] = Psi0

        def B(s, Psi, n):
            diag =  (Psi[n, 1:M-1 ] - Psi[n, 0:M-2 ])
            diag2 =  (Psi[n, 1:M-2 ] - Psi[n, 2:M-1 ])
            PsiMatrix = diags(np.array([ diag , diag2 ]), [0, 1], shape=(M-2, M-2)).toarray()
            if s == 1:
                return  -1 * a1 * PsiMatrix
            if s == -1:
                return a2 * PsiMatrix
        # Set up initial conditions for C

        rho1 = np.zeros([N, M])
        rho2 = np.zeros([N, M])
        Psi = np.zeros([N, M])
        E = np.zeros([N, M-1])

        rho1[0, :] = 0
        rho1[0, -1] = 1    

        rho2[0, :] = 0
        rho2[0, -1] = 1  

        Psi[0, :] = 0
        Psi[0, 0] = Psi0

        
        #Starting iteration
        for n in range(0, N-1):

             # Update border condition
            g1 = 1 / ( 1 + kf * dxi / ( D1 * kappa ) - (Psi[n,1]- Psi0))
            A1[0,0] = ( 1 - 2 * a1 + g1 * a1 )

            g2 = 1 / ( 1 + (Psi[n, 1] - Psi[n,0]))
            A2[0,0] = ( 1 - 2 * a2 + g2 * a2 )

            rho1[n+1, 1:M-1] = np.matmul(A1, rho1[n, 1:M-1])  + b1 + np.matmul(B(1, Psi, n), rho1[n, 1:M-1])
            rho1[n+1, 0] = g1 * rho1[n+1, 1]
            rho1[n+1, -1] = 1

            rho2[n+1, 1:M-1] = np.matmul(A2, rho2[n, 1:M-1]) + b2 + np.matmul(B(-1, Psi, n), rho2[n, 1:M-1]) 
            rho2[n+1, 0] = g2 * rho2[n+1, 1]
            rho2[n+1, -1] = 1

            Psi[n+1, 1:M-1] = np.matmul(Dinv, dxi * (rho2[n+1, 1:M-1] - rho1[n+1, 1:M-1]) -bPsi )
            Psi[n+1, 0] = Psi0
            Psi[n+1, -1] = 0
        
            E[n+1,0:M-1] = - self.kappa * (Psi[n+1,1:M] -  Psi[n+1,:M-1])/dxi
            
        print("Build Complete")
        self.rho1 = rho1
        self.rho2 = rho2
        self.Psi = Psi
        self.E = self.R * self.T * E/ (self.z * self.F)

        
    def remove_points(self, A, n):
        #n is the number of steps to skip
        if n >= 4:
            A = np.delete(A, [1, 2, 3])

        for i in range(0,int(len(A)/4)):
            index = i+n
            A = np.delete(A, [index-2, index-1, index])
        return A
    
    #Cm is the imported analytical solution
    def plot(self, tau, f, imageName='langmuir-diffusion-nernst'):
        
        self.imageName = imageName# + str(t)[-3:]
        Cb = self.Cb
        dtau = self.dtau
        C1 = Cb * self.rho1
        C2 = Cb * self.rho2

        if(f == 'E'):
            func = self.E
            ylabel = r'Electric Field (V/m)' 
            title = ''#'Electric Field In The Diffusion Problem With Nernst Interaction.'
        elif(f == 'phi'):
            func = self.R * self.T * self.Psi / (self.z * self.F)
            ylabel = r'Electric Potential (V)'
            title = ''#'Electric Potential In The Diffusion Problem With Nernst Interaction.'
             # this is done to avoid cluttering of numeric points over the analytic solution
        else:
            print("Unkown function. Need something real to plot")
            return -1
        
        kappa = self.kappa
        D1 = self.D1 
        mw = 4
        skip = 4
        
        xi2 = self.xi
        xi2 = xi2 / kappa * 1e9#* nanometerScale #change the scale of the scale to nanometer
        
        to_molar = 1e3
        for i in range(len(tau)):
            
            plt.figure(1)
            plt.style.use('thesis')

            fig, ax2 = plt.subplots()

            #color = 'tab:blue'
            ax2.tick_params(axis='y')#, labelcolor=color)
            ax1 = ax2.twinx() 
            #color = 'tab:red'

            plt.title(title, fontweight='bold')

            n = int(t[i]/dtau)
        
            ax1.plot(xi2, to_molar * C1[n, :], 'g^', label=r'$C_+$, $t =%.2f \mu s$' % Decimal(tau[i] * 1e6/ ( self.D1 * kappa ** 2 ) ))
            ax1.plot(xi2, to_molar * C2[n, :],'r^', label=r'$C_-$,  $t =%.2f \mu s$' % Decimal(tau[i] * 1e6/ ( self.D1 * kappa ** 2 ) ))
            ax1.legend(loc='upper left')
            if f == 'E':
                ax2.plot(xi2[0:self.M-1], func[n, :], 'b^' ,color='tab:blue', label=r'$\phi$,  $t = %.2f \mu s$' % Decimal(tau[i] * 1e6/ ( self.D1 * kappa ** 2 ) ))
            if f == 'phi':
                ax2.plot(xi2, func[n, :], 'b^', color='tab:blue', label=r'$\phi$,  $t =%.2f \mu s$' % Decimal(tau[i] * 1e6 / ( self.D1 * kappa ** 2 ) ))

            ax1.set_xlabel(r'Distance from the interface plate (nm)')
            ax1.set_ylabel(r'Molar Concentration')
            ax2.set_ylabel(ylabel)
            ax2.tick_params(axis='y')#, labelcolor=color)
            

            ax2.legend(loc=(0.01,  0.85))


            ################################## Plot parameters ##################################
            textstr = '\n'.join((
                r'$V_0=%.2f V$' % Decimal(self.V0),
                r'$k_f=%.2f m/s$' % Decimal(self.kf),
                r'$\kappa=%.2f nm^{-1}$' % Decimal(self.kappa /1e9),
                r'$C_b=%.2f M$' % Decimal(to_molar * self.Cb),))


            props = dict(boxstyle='round', facecolor='white', alpha=0.5)
            ax1.text(0.02, 0.8, textstr, fontsize=30, verticalalignment='top', bbox=props)
            
            
            fig.tight_layout()  # otherwise the right y-label is slightly clipped
            
            print("saving to:"+'../../img/'+ imageName + f + str(i) +'.eps')
            plt.savefig('../../img/'+ self.imageName + f + str(i) +'.eps')
    
            #plt.show()