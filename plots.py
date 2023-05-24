import numpy as np
import matplotlib.pyplot as plt
import funs

#####################################################
# Paper I: Dynamic view of the solid-state DNP effect
###################################################3
def plot_bloch(B1s,T1,T2):

    print("Paper I, Figure 5")

    R1 = 1/T1
    R2 = 1/T2


    fontsize = 16
    fig_scale = 1

    plt.rcParams['text.usetex'] = True


    nrows, ncols = 3, 4
    
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', sharey='row',
                           figsize=(fig_scale*(ncols*3.2+0.4),fig_scale*(nrows*2.0+0.4)))
    plt.subplots_adjust(wspace = 0.1,hspace=0.1)

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')


    for i,B1 in enumerate(B1s):
        print(f"B1:{B1} G")
       
        omega1 = funs.omega1_from_B1(B1)
        
        side = 30
        #scan offset frequency
        Omega = 2*np.pi * np.linspace(-65-side,140+side,2000) #M rad/s

        #EXACT
        fx, fy, fz = funs.bloch_ss(R1,R2,omega1,Omega)
        
        p =  R1 * fz
        sx =  omega1 * fx * p
        sy = -omega1 * fy * p

        for j in range(nrows):
            
            ax[j,i].axvline(x=0,linestyle=':',color='gray')
            ax[j,i].axvline(x=+14,linestyle='--',color='gray')
            ax[j,i].axvline(x=+45,linestyle='--',color='gray')
            ax[j,i].axvline(x=+140,linestyle='--',color='gray')
            
            ax[j,i].axhline(y=0,linestyle=':',color='gray')
            

            if j == 0:
                ax[j,i].plot(Omega/(2*np.pi),p,'-k')
                
                ax[j,i].set_title(f"$B_1=$ {B1} G",fontsize=fontsize) 
                
                ax[j,i].set_yticks([0,0.5,1])
                
                ax[j,i].axhline(y=1,linestyle=':',color='gray')
                
                if i ==0:
                    ax[j,i].set_ylabel(r"$s_z^{\rm ss}\, /\, s_z^{\rm eq}$",labelpad=20,fontsize=fontsize)
                
                
            elif j == 1:
                ax[j,i].plot(Omega/(2*np.pi),fx*omega1,'-',color='#ff7f0e',label="$x$")
                ax[j,i].plot(Omega/(2*np.pi),-fy*omega1,'-',color='#1f77b4',label="$y$")
                
                if i == 0:
                    ax[j,i].set_ylabel(r"$s_{x,y}^{\rm ss}\, /\, s_z^{\rm ss}$",labelpad=10,fontsize=fontsize)
                    ax[j,i].legend(loc = 'lower right',fontsize=fontsize-3, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)

                            
            elif j == 2:
                ax[j,i].plot(Omega/(2*np.pi),sx,'-',color='#ff7f0e',label="$x$")
                ax[j,i].plot(Omega/(2*np.pi),sy,'-',color='#1f77b4',label="$y$")
                
                ax[j,i].set_yticks([-0.15,0.0,0.15])
                                   
                ax[j,i].set_xticks([-50,0,14,45,100,140])
                ax[j,i].set_xticklabels([-50,0,"X","Q",100,"W"])#,fontsize=fontsize-2)
                if i == 0:
                    ax[j,i].set_ylabel(r"$s_{x,y}^{\rm ss}\, /\, s_z^{\rm eq}$",fontsize=fontsize)
                    ax[j,i].legend(loc = 'lower right',fontsize=fontsize-3, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)


    
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel("Offset frequency ($\Omega/2\\pi$)  [MHz]",fontsize=fontsize,labelpad=15)

    plt.savefig("pdfs/p1_fig05.pdf", bbox_inches='tight')
    return True



###################################################
def plot_Fxyz(B1s,freqIs,sides,T1,T2):

    print("Paper I, Figure 6b")

    R1 = 1/T1
    R2 = 1/T2

    fontsize = 16
    fig_scale = 1

    plt.rcParams['text.usetex'] = True

    nrows, ncols = 5, 4

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', sharey='row', 
                           figsize=(fig_scale*(ncols*3.2+0.4),fig_scale*(nrows*2.0+0.4)))
    plt.subplots_adjust(wspace = 0.1, hspace=0.1)


    for i,freqI in enumerate(freqIs):
        print(f"fI: {freqI} MHz")
        
        omega_I = 2*np.pi* freqI

        #scan offset frequency
        Omega = 2*np.pi * np.linspace(-(freqI+sides[i]),freqI+sides[i],2000) #M rad/s

        #EXACT
        B1 = B1s[-1]
        omega1 = funs.omega1_from_B1(B1)
        
        Fx, Fy, Fz, Fyp = funs.new_bloch_ss(R1,R2,omega1,omega_I,Omega)
        #approx
        Fz_app = funs.approxFz(R1,R2,omega1,omega_I,Omega)

        for j in range(nrows):
            
            ax[j,i].axvline(x=-freqI,linestyle='--',color='gray')
            ax[j,i].axvline(x=0,linestyle=':',color='gray')
            ax[j,i].axvline(x=+freqI,linestyle='--',color='gray')
            ax[j,i].axhline(y=0,linestyle=':',color='gray')

            ax[j,i].set_xticks([-freqI,0,freqI])
            
            if j == 0:
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.real(Fy),'-',color='#1f77b4',label="re")
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.imag(Fy),'--',color='#ff7f0e',label="im")
                
                if i == 0:
                    ax[j,i].set_ylabel("$F_y$  [ns]",fontsize=fontsize) 
                    ax[j,i].set_title("X band (9.2 GHz)",fontsize=fontsize) 
                if i == 1:
                    ax[j,i].set_title("K band (16.5 GHz)",fontsize=fontsize) 
                elif i == 2:
                    ax[j,i].set_title("Q band (30 GHz)",fontsize=fontsize) 
                elif i == 3:
                    ax[j,i].legend(loc = 'upper center',fontsize=fontsize-2, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)
                    ax[j,i].set_title("W band (92 GHz)",fontsize=fontsize) 
            
            elif j == 1:        
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.real(Fx),color='#1f77b4',label="re")
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.imag(Fx),'--',color='#ff7f0e',label="im")
                
                if i == 0:
                    ax[j,i].set_ylabel("$F_x$  [ns]",fontsize=fontsize) 
                
            elif j == 2:
                if i == 0:
                    ax[j,i].set_ylabel("$F_z$  [ns]",fontsize=fontsize) 
                    
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.real(Fz),color='#1f77b4',label="re")
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.imag(Fz),'--',color='#ff7f0e',label="im")
                
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.real(Fz_app),':k')#,label='approx')
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.imag(Fz_app),':k')
                
                if i == 3:
                    ax[j,i].legend(loc = 'upper right',fontsize=fontsize-2, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=2)
                  
            elif j == 3:
                
                B1 = B1s[-2] #G
                omega1 = funs.omega1_from_B1(B1)
                #EXACT
                Fx, Fy, Fz, Fyp = funs.new_bloch_ss(R1,R2,omega1,omega_I,Omega)
                #approx
                Fz_app = funs.approxFz(R1,R2,omega1,omega_I,Omega)
                
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.real(Fz),color='#1f77b4')#,label="re")
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.imag(Fz),'--',color='#ff7f0e')#,label="im")
                
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.real(Fz_app),':k',label='approx')
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.imag(Fz_app),':k')
                
                if i == 0:
                    ax[j,i].set_ylabel("$F_z$  [ns]",fontsize=fontsize)                
                elif i == 3:
                    ax[j,i].legend(loc = 'upper right',fontsize=fontsize-2, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)
                
                    
            elif j == 4:
                
                B1 = B1s[-3] #G
                omega1 = funs.omega1_from_B1(B1)
                #EXACT
                Fx, Fy, Fz, Fyp = funs.new_bloch_ss(R1,R2,omega1,omega_I,Omega)
                #approx
                Fz_app = funs.approxFz(R1,R2,omega1,omega_I,Omega)
                
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.real(Fz),color='#1f77b4',label="re")
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.imag(Fz),'--',color='#ff7f0e',label="im")
                
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.real(Fz_app),':k',label='approx')
                ax[j,i].plot(Omega/(2*np.pi),1e3*np.imag(Fz_app),':k')
                if i == 0:
                        ax[j,i].set_ylabel("$F_z$  [ns]",fontsize=fontsize)
            if j > 1:
                ax[j,i].set_ylim(-22,22)
                
                
    fig.text(0.05, 0.50,f'$B_1={B1s[-1]}$ G', va='center', rotation='vertical',fontsize=fontsize)
    fig.text(0.05, 0.345,f'$B_1={B1s[-2]}$ G', va='center', rotation='vertical',fontsize=fontsize)
    fig.text(0.05, 0.19,f'$B_1={B1s[-3]}$ G', va='center', rotation='vertical',fontsize=fontsize)

    
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel("Offset frequency ($\Omega/2\\pi$)  [MHz]",fontsize=fontsize,labelpad=15)

    plt.savefig("pdfs/p1_fig06b.pdf", bbox_inches='tight')
    return True



#########################
def plot_newbloch(B1,freqIs,sides,T1,T2):

    print("Paper I, Figure 7")

    omega1 = funs.omega1_from_B1(B1)
    R1 = 1/T1
    R2 = 1/T2

    fontsize = 16
    fig_scale = 1

    plt.rcParams['text.usetex'] = True

    nrows, ncols = 4, 4

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', 
                           figsize=(fig_scale*(ncols*3.2+0.4),fig_scale*(nrows*2.0+0.4)))
    plt.subplots_adjust(wspace = 0.2)


    for i,freqI in enumerate(freqIs):
        print(f"fI: {freqI} MHz")
        
        omega_I = 2*np.pi* freqI

        #scan offset frequency
        Omega = 2*np.pi * np.linspace(-(freqI+sides[i]),freqI+sides[i],2000) #M rad/s

        #EXACT
        Ti, Ti0, vp_d2, Tz, Txp, Typ, Tx = funs.TiTz_exact(R1,R2,omega1,omega_I,Omega)
        Ti0 = Ti - vp_d2
        #Bloch
        fx, fy, fz = funs.bloch_ss(R1,R2,omega1,Omega)
        
        sz =  R1 * fz
        sx_sz =  omega1 * fx 
        sy_sz = -omega1 * fy 

        TX = omega1*fx*Txp
        TY = -omega1*fy*Typ
        #Tis = TisX + TisY
        
        for j in range(nrows):
            
            ax[j,i].axvline(x=-freqI,linestyle='--',color='gray')
            ax[j,i].axvline(x=0,linestyle=':',color='gray')
            ax[j,i].axvline(x=+freqI,linestyle='--',color='gray')
            ax[j,i].axhline(y=0,linestyle=':',color='gray')
            
            ax[j,i].set_xticks([-freqI,0,freqI]) 
            
            if j == 0:
                ax[j,i].plot(Omega/(2*np.pi),sx_sz,'-',color='C1',label="$x$")
                ax[j,i].plot(Omega/(2*np.pi),sy_sz,'-',color='C0',label="$y$")
                
                if i == 0:
                    ax[j,i].set_ylabel(r"$s_{x,y}^{\rm ss}\, /\, s_z^{\rm ss}$",labelpad=10,fontsize=fontsize)
                    ax[j,i].set_title("X band (9.2 GHz)",fontsize=fontsize)   
                if i == 1:
                    ax[j,i].set_title("K band (16.5 GHz)",fontsize=fontsize) 
                elif i == 2:
                    ax[j,i].set_title("Q band (30 GHz)",fontsize=fontsize) 
                elif i == 3:
                    ax[j,i].legend(loc = 'lower right',fontsize=fontsize-2, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)
                    ax[j,i].set_title("W band (92 GHz)",fontsize=fontsize) 
                elif i == 4:
                    ax[j,i].set_title("J band",fontsize=fontsize) 
            
            elif j == 1:        
                ax[j,i].plot(Omega/(2*np.pi),1e3*Txp,color='C1',label="$x$")
                ax[j,i].plot(Omega/(2*np.pi),1e3*Typ,color='C0',label="$y$")
                
                if i == 0:
                    ax[j,i].set_ylabel(r"$T_{x,y}'$  [ns]",fontsize=fontsize)            
                if i == 3:
                    ax[j,i].legend(loc = 'upper center',fontsize=fontsize-3, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)
                    
            elif j == 2:
                ax[j,i].plot(Omega/(2*np.pi),1e3*Ti,'-r')
                ax[j,i].plot(Omega/(2*np.pi),1e3*Ti0,'--r',label=r"$T_i^0$")
                if i == 0:
                    ax[j,i].set_ylabel(r"$T_i$  [ns]",labelpad=10,fontsize=fontsize)
                if i == 3:
                    ax[j,i].legend(loc = 'upper center',fontsize=fontsize-3, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)
                    
            elif j == 3:            
                ax[j,i].plot(Omega/(2*np.pi),1e3*TX,color='C1',label="$x$")
                ax[j,i].plot(Omega/(2*np.pi),1e3*TY,color='C0',label="$y$")
                
                ax[j,i].plot(Omega/(2*np.pi),1e3*Tz,'--',color='black')
                if i == 0:
                    ax[j,i].set_ylabel(r"$T_z$  [ns]",labelpad=0,fontsize=fontsize)
                if i == 3:
                    ax[j,i].legend(loc = 'upper center',fontsize=fontsize-3, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)
                   
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel("Offset frequency ($\Omega/2\\pi$)  [MHz]",fontsize=fontsize,labelpad=15)

    plt.savefig(f"pdfs/p1_fig07_B{B1}.pdf", bbox_inches='tight')

    #plt.show()
    return True


##############################
def plot_rates(B1,freqIs,sides,T1,T2):

    print("Paper I, Figure 8")

    omega1 = funs.omega1_from_B1(B1)
    R1 = 1/T1
    R2 = 1/T2

    fontsize = 16
    fig_scale = 1

    plt.rcParams['text.usetex'] = True

    nrows, ncols = 3, 4

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', 
                           figsize=(fig_scale*(ncols*3.2+0.4),fig_scale*(nrows*2.0+0.4)))
    plt.subplots_adjust(wspace = 0.2)

    for i,freqI in enumerate(freqIs):

        print(f"fI: {freqI} MHz")
        
        omega_I = 2*np.pi* freqI

        #scan offset frequency
        Omega = 2*np.pi * np.linspace(-(freqI+sides[i]),freqI+sides[i],2000) #M rad/s

        #EXACT
        Ti, Ti0, vp_d2, Tz, Txp, Typ, Tx = funs.TiTz_exact(R1,R2,omega1,omega_I,Omega)
        #Classical
        vZ, vD = funs.vZvD_classical(R1,R2,omega1,omega_I,Omega)

        #Bloch
        fx, fy, p = funs.bloch_ss(R1,R2,omega1,Omega)
        sz =  R1 * p
        
        for j in range(nrows):
            
            ax[j,i].axvline(x=-freqI,linestyle='--',color='gray')
            ax[j,i].axvline(x=0,linestyle=':',color='gray')
            ax[j,i].axvline(x=+freqI,linestyle='--',color='gray')
            ax[j,i].axhline(y=0,linestyle=':',color='gray')
            
            ax[j,i].set_xticks([-freqI,0,freqI]) 
            
            if j == 0:
                ax[j,i].plot(Omega/(2*np.pi),1e3*vp_d2,'-r',label='exact')#,color='black')
                ax[j,i].plot(Omega/(2*np.pi),1e3*(vZ+vD),'--k',label='classical')
                
                if i == 0:
                    ax[j,i].set_ylabel(r"$v_+/\delta^2$  [ns]",labelpad=10,fontsize=fontsize)
                    ax[j,i].set_title("X band (9.2 GHz)",fontsize=fontsize) 
                    
                if i == 1:
                    ax[j,i].set_title("K band (16.5 GHz)",fontsize=fontsize) 
                elif i == 2:
                    ax[j,i].set_title("Q band (30 GHz)",fontsize=fontsize) 
                elif i == 3:
                    ax[j,i].legend(loc = 'upper center',fontsize=fontsize-3, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)
                    ax[j,i].set_title("W band (92 GHz)",fontsize=fontsize) 
                elif i == 4:
                    ax[j,i].set_title("J band",fontsize=fontsize) 
            
            elif j == 1:
                ax[j,i].plot(Omega/(2*np.pi),1e3*Tz,'-',color='#2ca02c',label='exact')
                ax[j,i].plot(Omega/(2*np.pi),1e3*(vD-vZ),'--k',label='classical')

                
                if i == 0:
                    ax[j,i].set_ylabel(r"$v_-/\delta^2$  [ns]",fontsize=fontsize)            
                if i == 3:
                    ax[j,i].legend(loc = 'upper center',fontsize=fontsize-3, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)
                    
            elif j == 2:
                ax[j,i].plot(Omega/(2*np.pi),1e3*(Tz+vp_d2)/2,'-',color='#8c564b',label='exact')
                ax[j,i].plot(Omega/(2*np.pi),1e3*vD,'--k',label='classical')

                if i == 0:
                    ax[j,i].set_ylabel(r"$v_2/\delta^2$  [ns]",labelpad=8,fontsize=fontsize)
                if i == 3:
                    ax[j,i].legend(loc = 'upper center',fontsize=fontsize-3, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)
      
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel("Offset frequency ($\Omega/2\\pi$)  [MHz]",fontsize=fontsize,labelpad=15)


    plt.savefig(f"pdfs/p1_fig08_B{B1}.pdf", bbox_inches='tight')
    return True





#####################
def plot_enhancement(B1,freqIs,sides,T1,T2,T1n,closest_distance_nm,radical_conc_molar):

    print("Paper I, Figure 9")

    omega1 = funs.omega1_from_B1(B1)
    R1 = 1/T1
    R2 = 1/T2

    R1n = 1/T1n
    delta2 = funs.delta2(closest_distance_nm, radical_conc_molar)

    R1n_del2 = R1n/delta2
    print(1e3*R1n_del2," ns")

    fontsize = 16
    fig_scale = 1

    plt.rcParams['text.usetex'] = True


    nrows, ncols = 3, 4

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', #sharey='row',
                           figsize=(fig_scale*(ncols*3.2+0.4),fig_scale*(nrows*2.0+0.4)))
    plt.subplots_adjust(wspace = 0.2)#,hspace=0.1)

    for i,freqI in enumerate(freqIs):

        print(f"fI: {freqI} MHz")
        
        omega_I = 2*np.pi* freqI

        #scan offset frequency
        Omega = 2*np.pi * np.linspace(-(freqI+sides[i]),freqI+sides[i],2000) #M rad/s

        #Bloch
        fx, fy, fz = funs.bloch_ss(R1,R2,omega1,Omega)
        #EXACT
        Ti, Ti0, vp_d2, Tz, Txp, Typ, Tx  = funs.TiTz_exact(R1,R2,omega1,omega_I,Omega)
        #Classical
        vZ, vD = funs.vZvD_classical(R1,R2,omega1,omega_I,Omega)
        
        p =  R1 * fz
        
        pvm_Rdel = p*Tz/R1n_del2
        pX = R1n_del2/(R1n_del2 + vp_d2)
        epsSE = pvm_Rdel*pX*658
        
        pvm_Rdel_cl = p*(vD-vZ)/R1n_del2
        pX_cl = R1n_del2/(R1n_del2 + (vD+vZ))   
        epsSE_cl = pvm_Rdel_cl*pX_cl*658


        for j in range(nrows):
            
            ax[j,i].axvline(x=-freqI,linestyle='--',color='gray')
            ax[j,i].axvline(x=0,linestyle=':',color='gray')
            ax[j,i].axvline(x=+freqI,linestyle='--',color='gray')
            ax[j,i].axhline(y=0,linestyle=':',color='gray')

            ax[j,i].set_xticks([-freqI,0,freqI]) 
 
            if j == 0:
                ax[j,i].plot(Omega/(2*np.pi),pX,'-',color='#d62728',label="exact")
                ax[j,i].plot(Omega/(2*np.pi),pX_cl,'--k',label="classical")
                
                ax[j,i].axhline(y=1,linestyle=':',color='gray')
                ax[j,i].set_yticks([0,0.5,1])
                ax[j,i].set_ylim(-0.05,1.05)
                
                if i == 0:
                    ax[j,i].set_title("X band (9.2 GHz)",fontsize=fontsize) 
                    ax[j,i].set_ylabel(r"$p_X$",fontsize=fontsize,labelpad=10)
                if i == 1:
                    ax[j,i].set_title("K band (16.5 GHz)",fontsize=fontsize) 
                elif i == 2:
                    ax[j,i].set_title("Q band (30 GHz)",fontsize=fontsize) 
                elif i == 3:
                    ax[j,i].legend(loc = 'lower right',fontsize=fontsize-3, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)
                    ax[j,i].set_title("W band (92 GHz)",fontsize=fontsize) 
                elif i == 4:
                    ax[j,i].set_title("J band",fontsize=fontsize) 
            
            elif j == 1:
                ax[j,i].plot(Omega/(2*np.pi),pvm_Rdel,color='#2ca02c',label="exact")
                ax[j,i].plot(Omega/(2*np.pi),pvm_Rdel_cl,'--k',label="classical")#,color='gray')
                        
                if i == 0:
                    ax[j,i].set_ylabel(r"$pv_-/R_{1I}$",fontsize=fontsize,labelpad=10)
                if i == 3:
                    ax[j,i].legend(loc = 'lower right',fontsize=fontsize-3, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)
                
            elif j == 2:
                ax[j,i].plot(Omega/(2*np.pi),epsSE,'-',color='#9467bd',label="exact") #8c564b #bcbd22
                ax[j,i].plot(Omega/(2*np.pi),epsSE_cl,'--k',label="classical")
                
                if i == 0:
                    ax[j,i].set_ylabel(r"$\epsilon_{SE}$",fontsize=fontsize)
                if i == 3:
                    ax[j,i].legend(loc = 'lower right',fontsize=fontsize-3, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)
                
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel("Offset frequency ($\Omega/2\\pi$)  [MHz]",fontsize=fontsize,labelpad=15)

    plt.savefig(f"pdfs/p1_fig09_B{B1}.pdf", bbox_inches='tight')
    return True


##############################
def plot_filters(B1,freqIs,sides,T1,T2,choice):

    omega1 = funs.omega1_from_B1(B1)
    R1 = 1/T1
    R2 = 1/T2

    fontsize = 14
    fig_scale = 1

    plt.rcParams['text.usetex'] = True

    #calculation
    freqI = 400
    sides = 300

    freqI = 140
    sides = 60

    omega_I = 2*np.pi* freqI


    fig = plt.figure(figsize=(4,5))

    #scan offset frequency
    pntN = 431
    Omega = 2*np.pi * np.linspace(-freqI-sides,freqI+sides,pntN) #M rad/s
    
    #Bloch
    omega1 = funs.omega1_from_B1(B1)
    fx, fy, p = funs.bloch_ss(R1,R2,omega1,Omega)
    sz =  R1 * p
    sx =  omega1 * fx * sz
    sy =  omega1 * fy * sz

    Fx, Fy, Fz, Fyp = funs.new_bloch_ss(R1,R2,omega1,omega_I,Omega)

    Tx = np.real(1j*Fz*omega1*Fyp)


    if choice == 'class':
        print("Paper I, Figure 10a")
        sz_norm = sz/np.max(sz)
        v = omega1*fx*Tx
        v_norm = (v)/np.max(v)
        
        plt.plot(Omega/(2*np.pi),sz+4,'-',color='k',linewidth=2)
        plt.plot(Omega/(2*np.pi),1500*v+1,'-',color='C2',linewidth=2)
        #plt.plot(Omega/(2*np.pi),prod_norm,color='C5',label="DNP enhancement",linewidth=2)

    else:
        print("Paper I, Figure 10b")
        sx_norm = sx/np.max(sx)
        Tx_norm = Tx/np.max(Tx)
        #prod = sx_norm*Tx_norm
        #prod_norm = prod/np.max(prod)

        plt.plot(Omega/(2*np.pi),0.4*sx_norm+4,'-',color='C0',linewidth=2)#,color='#ff7f0e'
        plt.plot(Omega/(2*np.pi),600*Tx+1,'-',color='C1',linewidth=2)
        #plt.plot(Omega/(2*np.pi),prod_norm,color='C5',label="DNP enhancement",linewidth=2)


    plt.ylim(0,6)

    plt.xlabel('Offset frequency',labelpad=10,fontsize=fontsize-2)

    label_p = r"$+\omega_I$"
    label_m = r"$-\omega_I$"
    plt.xticks([-freqI,0,freqI],[label_m,0,label_p],fontsize=fontsize-2)

    plt.yticks([])
    #plt.axis('off')
    for pos in ['right', 'top', 'bottom', 'left']:
        plt.gca().spines[pos].set_visible(False)

    if choice == 'class':
        plt.savefig(f"pdfs/p1_fig10a.pdf", bbox_inches='tight')

    else:
        plt.savefig(f"pdfs/p1_fig10b.pdf", bbox_inches='tight')

    return True

###################################################
# Paper II: The solid-state DNP effect in liquids
##################################################
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#########################
def plot_liquid(plot,B1,freqIs,sides,T1,T2,tau,delta2,T1n_min):

    omega1 = funs.omega1_from_B1(B1)
    R1 = 1/T1
    R2 = 1/T2

    fontsize = 16
    fig_scale = 1

    plt.rcParams['text.usetex'] = True

    nrows, ncols = 6, len(freqIs)

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', 
                           figsize=(fig_scale*(ncols*3.2+0.4),fig_scale*(nrows*2.0+0.4)))
    plt.subplots_adjust(wspace = 0.2)

    for i,freqI in enumerate(freqIs):
        
        omega_I = 2*np.pi* freqI

        #scan offset frequency
        Omega = 2*np.pi * np.linspace(-(freqI+sides[i]),freqI+sides[i],2000) #M rad/s

        #SOLID
        Tii_s,Tii0_s,vp_s,Tiz_s,Tix_s,Tiy_s,Tixp_s = funs.TiTz_exact(R1,R2,omega1,omega_I,Omega)
        #Numerical 
        Tii_l,Tii0_l,vp_l, Tix_l, Tiy_l = funs.eig_numerical(R1,R2,omega1,omega_I,Omega,'exp',tau)
        Tii_f,Tii0_f,vp_f, Tix_f, Tiy_f = funs.eig_numerical(R1,R2,omega1,omega_I,Omega,'ffhs',tau)
        #APPROX 
        Tii_p, Tii0_p, vp_p, Tixp_p = funs.eig_perturbative(R1,R2,omega1,omega_I,Omega,'ffhs',tau)
        #Classical
        #Tii_cl, Tis_cl = TiiTis_classical(R1,R2,omega1,omega_I,Omega)

        R1n = np.max([2*delta2*Tii0_l,1/T1n_min])
        print(1e-3/R1n," ms")
        R1n_del2 = R1n/delta2
        print(1e3*R1n_del2," ns")

        print("freq:",freqI)
        print(R1n_del2-Tii0_s)
        print(R1n_del2-Tii0_l)
        print(R1n_del2-Tii0_f)
        print(R1n_del2-Tii0_p,"\n")
        
        #Bloch
        fx, fy, p = funs.bloch_ss(R1,R2,omega1,Omega)
        
        sz =  R1 * p
        sx =  omega1 * fx * sz
        sy = -omega1 * fy * sz

        TisX = omega1*fx*Tix_l
        TisY = -omega1*fy*Tiy_l
        Tiz_l = TisX + TisY
        Tixp_l = Tiz_l/omega1/fx

        TisX = omega1*fx*Tix_f
        TisY = -omega1*fy*Tiy_f
        Tiz_f = TisX + TisY
        Tixp_f = Tiz_f/omega1/fx

        Tiz_p = omega1*fx*Tixp_p
        
        for j in range(nrows):
            
            ax[j,i].axvline(x=-freqI,linestyle='--',color='gray')
            ax[j,i].axvline(x=0,linestyle=':',color='gray')
            ax[j,i].axvline(x=+freqI,linestyle='--',color='gray')

            if j<4:
                ax[j,i].axhline(y=0,linestyle=':',color='gray')
            
            ax[j,i].set_xticks([-freqI,0,freqI]) 
            
            if j == 0:
                ax[j,i].plot(Omega/(2*np.pi),sx,'-',color='#1f77b4',label="$x$")

                if i == 0:
                    ax[j,i].set_ylabel(r"$s_{x}^{\rm ss}\, /\, s_z^{\rm eq}$",fontsize=fontsize)
                    ax[j,i].set_title("X band (9.2 GHz)",fontsize=fontsize) 
                elif i == 5:
                    ax[j,i].set_title("K band (16.5 GHz)",fontsize=fontsize) 
                elif i == 1:
                    ax[j,i].set_title("Q band (30 GHz)",fontsize=fontsize) 
                elif i == 2:
                    ax[j,i].set_title("W band (92 GHz)",fontsize=fontsize) 
                elif i == 3:
                    ax[j,i].set_title("J band (260 GHz)",fontsize=fontsize) 
            
            elif j == 1:  
                if plot == 'solid':
                    ax[j,i].plot(Omega/(2*np.pi),1e3*Tixp_s,color='#ff7f0e',label="solid") 
                    ax[j,i].plot(Omega/(2*np.pi),1e3*Tixp_l,'--k',label="liquid")
                elif plot == 'liquid':
                    ax[j,i].plot(Omega/(2*np.pi),1e6*Tixp_l,'--k',label="exp")
                    ax[j,i].plot(Omega/(2*np.pi),1e6*Tixp_f,color='#ff7f0e',label="ffhs") 
                else:
                    ax[j,i].plot(Omega/(2*np.pi),1e6*Tixp_f,color='#ff7f0e',label="ffhs")
                    ax[j,i].plot(Omega/(2*np.pi),1e6*Tixp_p,'-.k',label="pert") 
                
                if i == 0:
                    if plot == 'solid':
                        ax[j,i].set_ylabel(r"$T_x$  [ns]",labelpad=10,fontsize=fontsize)  
                    else:
                        ax[j,i].set_ylabel(r"$T_x$  [ps]",labelpad=10,fontsize=fontsize)  
                if i == 3:
                    ax[j,i].legend(loc = 'upper center',fontsize=fontsize-3, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)
                    
            elif j == 2:

                if plot == 'solid':
                    ax[j,i].plot(Omega/(2*np.pi),1e3*vp_s,'-r')#,label='solid')
                    ax[j,i].plot(Omega/(2*np.pi),1e3*vp_l,'--k')#,label='liquid')

                    ax[j,i].axhline(y=1e3*Tii0_s,linestyle=':',color='r',label=r"$T_i^0$ solid")
                    ax[j,i].axhline(y=1e3*Tii0_l,linestyle=':',color='k',label=r"$T_i^0$ liquid")
                elif plot == 'liquid':
                    ax[j,i].plot(Omega/(2*np.pi),1e6*vp_l,'--k',label='exp')
                    ax[j,i].plot(Omega/(2*np.pi),1e6*vp_f,'-r',label='ffhs')
                else:
                    ax[j,i].plot(Omega/(2*np.pi),1e6*vp_f,'-r',label='ffhs')
                    ax[j,i].plot(Omega/(2*np.pi),1e6*vp_p,'-.k',label='pert')
                if i == 0:
                    if plot == 'solid':
                        ax[j,i].set_ylabel(r"$v_+/\langle\delta^2\rangle$  [ns]",labelpad=10,fontsize=fontsize)
                    else:
                        ax[j,i].set_ylabel(r"$v_+/\langle\delta^2\rangle$  [ps]",labelpad=10,fontsize=fontsize)
                if plot == 'pert':
                    if i == 3:
                        ax[j,i].legend(loc = 'upper center',fontsize=fontsize-3, 
                                  facecolor='white',framealpha=1,frameon=True,ncol=1)
                else:
                    if i == 3:
                        ax[j,i].legend(loc = 'upper center',fontsize=fontsize-3, 
                                  facecolor='white',framealpha=1,frameon=True,ncol=1)
                 
            elif j == 3:

                if plot == 'solid':
                    ax[j,i].plot(Omega/(2*np.pi),sz*1e3*Tiz_s,color='#2ca02c',label='solid')
                    ax[j,i].plot(Omega/(2*np.pi),sz*1e3*Tiz_l,'--k',label='liquid')
                elif plot == 'liquid':
                    ax[j,i].plot(Omega/(2*np.pi),sz*1e6*Tiz_l,'--k',label='exp')
                    ax[j,i].plot(Omega/(2*np.pi),sz*1e6*Tiz_f,color='#2ca02c',label='ffhs')
                else:
                    ax[j,i].plot(Omega/(2*np.pi),sz*1e6*Tiz_f,color='#2ca02c',label='ffhs')
                    ax[j,i].plot(Omega/(2*np.pi),sz*1e6*Tiz_p,'-.k',label='pert')
                if i == 0:
                    if plot == 'solid':
                        ax[j,i].set_ylabel(r"$p\,v_-/\langle\delta^2\rangle$  [ns]",fontsize=fontsize)
                    else:
                        ax[j,i].set_ylabel(r"$p\,v_-/\langle\delta^2\rangle$  [ps]",fontsize=fontsize)
                if i == 3:
                    ax[j,i].legend(loc = 'lower right',fontsize=fontsize-3, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)   

            elif j == 4:
                if plot == 'solid':
                    ax[j,i].plot(Omega/(2*np.pi),658*sz*Tiz_s/(R1n_del2+vp_s),color='#9467bd',label='solid')
                    ax[j,i].plot(Omega/(2*np.pi),658*sz*Tiz_l/(R1n_del2+vp_l),'--k',label='liquid')
                elif plot == 'liquid':
                    ax[j,i].plot(Omega/(2*np.pi),658*sz*Tiz_l/(R1n_del2+vp_l),'--k',label='exp')
                    ax[j,i].plot(Omega/(2*np.pi),658*sz*Tiz_f/(R1n_del2+vp_f),color='#9467bd',label='ffhs')
                else:
                    ax[j,i].plot(Omega/(2*np.pi),658*sz*Tiz_f/(R1n_del2+vp_f),color='#9467bd',label='ffhs')
                    ax[j,i].plot(Omega/(2*np.pi),658*sz*Tiz_p/(R1n_del2+vp_p),'-.k',label='pert')

                if i == 0:
                    ax[j,i].set_ylabel(r"$\epsilon_{\rm SE}$",labelpad=10,fontsize=fontsize)
                if i == 3:
                    ax[j,i].legend(loc = 'lower right',fontsize=fontsize-3, 
                                   facecolor='white',framealpha=1,frameon=True,ncol=1)
                          
            elif j == 5:

                if plot == 'solid':
                    ax[j,i].plot(Omega/(2*np.pi),R1n_del2/(R1n_del2+vp_s),color='#d62728',label='solid')
                    ax[j,i].plot(Omega/(2*np.pi),R1n_del2/(R1n_del2+vp_l),'--k',label='liquid')
                elif plot == 'liquid':
                    ax[j,i].plot(Omega/(2*np.pi),R1n_del2/(R1n_del2+vp_l),'--k',label='exp')
                    ax[j,i].plot(Omega/(2*np.pi),R1n_del2/(R1n_del2+vp_f),color='#d62728',label='ffhs')
                else:
                    ax[j,i].plot(Omega/(2*np.pi),R1n_del2/(R1n_del2+vp_f),color='#d62728',label='ffhs')
                    ax[j,i].plot(Omega/(2*np.pi),R1n_del2/(R1n_del2+vp_p),'-.k',label='pert')
                
                ax[j,i].axhline(y=1,linestyle=':',color='gray')
                if i == 0:
                    ax[j,i].set_ylabel(r"$p_X$",labelpad=15,fontsize=fontsize)
                if plot == 'pert':
                    if i == 3:
                        ax[j,i].legend(loc = 'lower center',fontsize=fontsize-3, 
                                  facecolor='white',framealpha=1,frameon=True,ncol=1)
                else:
                    if i == 3:
                        ax[j,i].legend(loc = 'lower center',fontsize=fontsize-3, 
                                  facecolor='white',framealpha=1,frameon=True,ncol=1)
            
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel("Offset frequency ($\Omega/2\\pi$)  [MHz]",fontsize=fontsize,labelpad=15)

    if plot == 'solid':
        print("Paper II, Figure 3")
        plt.savefig(f"pdfs/p2_fig03_B{B1}.pdf", bbox_inches='tight')
    elif plot == 'liquid':
        print("Paper II, Figure 4")
        plt.savefig(f"pdfs/p2_fig04_B{B1}.pdf", bbox_inches='tight')
    elif plot == 'pert':
        print("Paper II, Figure 5")
        plt.savefig(f"pdfs/p2_fig05_B{B1}.pdf", bbox_inches='tight')
    return True

##############
def plot_Tx(freqI,sides,B1,T1e,T2e,tau):

    print("Paper II, Figure 6")

    omega_I = 2*np.pi* freqI
    omega1 = funs.omega1_from_B1(B1)
    R1 = 1/T1e
    R2 = 1/T2e


    #1) calculate Tix' at Omega = omega_I 
    Omega = np.array([omega_I])

    #first the solid baseline
    #SOLID
    Tii_s,Tii0_s,vp_s,Tiz_s,Tix_s,Tiy_s,Tixp_s = funs.TiTz_exact(R1,R2,omega1,omega_I,Omega)
    solid0 = Tixp_s

    #scan motional time scale relative to T2e
    tau_T2_0 = np.logspace(-2,4,200)
    
    exp0 = np.zeros_like(tau_T2_0)
    ffhs0 = np.zeros_like(tau_T2_0)

    for i,ratio in enumerate(tau_T2_0):
        tau = ratio/R2

        #Numerical 
        Tii_e,Tii0_e,vp_e, Tix_e, Tiy_e = funs.eig_numerical(R1,R2,omega1,omega_I,Omega,'exp',tau)
        Tii_f,Tii0_f,vp_f, Tix_f, Tiy_f = funs.eig_numerical(R1,R2,omega1,omega_I,Omega,'ffhs',tau)
        
        fx, fy, p = funs.bloch_ss(R1,R2,omega1,Omega)
        sz =  R1 * p
        sx =  omega1 * fx * sz
        sy = -omega1 * fy * sz

        TisX = omega1*fx*Tix_e
        TisY = -omega1*fy*Tiy_e
        Tiz = TisX + TisY
        Tixp_e = Tiz/omega1/fx

        TisX = omega1*fx*Tix_f
        TisY = -omega1*fy*Tiy_f
        Tiz = TisX + TisY
        Tixp_f = Tiz/omega1/fx
        
        exp0[i] = Tixp_e/solid0
        ffhs0[i] = Tixp_f/solid0


    #2) calculate Tix' at Omega = omega1*np.sqrt(R2/R1)
    Omega = np.array([omega1*np.sqrt(R2/R1)])

    #first the solid baseline
    #SOLID
    Tii_s,Tii0_s,vp_s,Tiz_s,Tix_s,Tiy_s,Tixp_s = funs.TiTz_exact(R1,R2,omega1,omega_I,Omega)
    solid1 = Tixp_s

    #scan motional time scale relative to T2e
    tau_T2_1 = np.logspace(-3,2,200)
    
    exp1 = np.zeros_like(tau_T2_1)
    ffhs1 = np.zeros_like(tau_T2_1)
    
    for i,ratio in enumerate(tau_T2_1):

        tau = ratio/R2
        #Numerical 
        Tii_e,Tii0_e,vp_e, Tix_e, Tiy_e = funs.eig_numerical(R1,R2,omega1,omega_I,Omega,'exp',tau)
        Tii_f,Tii0_f,vp_f, Tix_f, Tiy_f = funs.eig_numerical(R1,R2,omega1,omega_I,Omega,'ffhs',tau)
        
        fx, fy, fz = funs.bloch_ss(R1,R2,omega1,Omega)
        p =  R1 * fz
        sx =  omega1 * fx * p
        sy = -omega1 * fy * p

        TisX = omega1*fx*Tix_e
        TisY = -omega1*fy*Tiy_e
        Tiz = TisX + TisY
        Tixp_e = Tiz/omega1/fx

        TisX = omega1*fx*Tix_f
        TisY = -omega1*fy*Tiy_f
        Tiz = TisX + TisY
        Tixp_f = Tiz/omega1/fx
        
        exp1[i] = Tixp_e/solid1
        ffhs1[i] = Tixp_f/solid1

    #print(solid)
    fontsize = 14
    fig_scale = 1

    plt.rcParams['text.usetex'] = True

    nrows, ncols = 2, 1

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, #sharex='col', #sharey='row',
                           figsize=(fig_scale*(ncols*3.8+0.4),fig_scale*(nrows*2.8+0.4)))
    plt.subplots_adjust(hspace = 0.2)#,hspace=0.1)

    for j in range(nrows):

        if j == 0:
            ax[j].set_xscale("log")
            ax[j].axvline(x=1,linestyle=':',color='gray')
            ax[j].axhline(y=0,linestyle=':',color='gray')
            ax[j].hlines(y=0.5,xmin=0,xmax=10,linestyle=':',color='gray')
            ax[j].axhline(y=1,linestyle=':',color='gray')

            ax[j].plot(tau_T2_0,exp0,'--',label="exp",color='C1')
            ax[j].plot(tau_T2_0,ffhs0,'-',label="ffhs",color="C0")
            
            ax[j].set_xticks([1e-2,1e-1,1e0,1e1,1e2,1e3,1e4])

            ax[j].set_ylabel("Reduction factor",fontsize=fontsize)

            ax[j].legend(frameon=False,loc="upper left",fontsize=fontsize-2)

            inset_ax = inset_axes(ax[j],
                    width="35%", # width = 30% of parent_bbox
                    height="40%", # height : 1 inch
                    loc="right")
            inset_ax.set_xlim(-0.5,8)
            inset_ax.set_xticks([0,4,8])
            inset_ax.axvline(x=1,linestyle=':',color='gray')
            inset_ax.axhline(y=0,linestyle=':',color='gray')
            inset_ax.axhline(y=0.5,linestyle=':',color='gray')
            inset_ax.axhline(y=1,linestyle=':',color='gray')
            inset_ax.plot(tau_T2_0,ffhs0,'-',color="C0")
            inset_ax.plot(tau_T2_0,exp0,'--',color='C1')

        elif j == 1:
            ax[j].set_xscale("log")
            ax[j].set_yscale("log")
            ax[j].axvline(x=1,linestyle=':',color='gray')
            ax[j].axhline(y=1,linestyle=':',color='gray')

            ax[j].plot(tau_T2_1,exp1,'--',label="exp",color='C1')
            ax[j].plot(tau_T2_1,ffhs1,'-',label="ffhs",color="C0")
            
            ax[j].set_xticks([1e-3,1e-2,1e-1,1e0,1e1,1e2])
            #ax[j].set_yticks([1e-1,1e0,1e1,1e2,1e3])
            ax[j].minorticks_off()

            ax[j].set_xlabel(r"Relative motional time scale $\tau/T_{2S}$",fontsize=fontsize)
            ax[j].set_ylabel("Magnification factor",fontsize=fontsize)

            ax[j].legend(frameon=False,loc="lower left",fontsize=fontsize-2)

            
    plt.savefig(f"pdfs/p2_fig06.pdf", bbox_inches='tight')
    return True

#################################
# Comparison with experiment
###############################

def plot_epr(epr_J,T2e,phase,scale):

    print("Paper 2, cw-EPR\n")

    freq_axis = funs.field_to_frequency(epr_J[:,0],'epr')
    Lorentz_deriv_re, Lorentz_deriv_im  = funs.Lorentzian_derivative(freq_axis,0,1/T2e)

    phase_rad = np.deg2rad(phase)
    cos = np.cos(phase_rad)
    sin = np.sin(phase_rad)
    Lorentz_deriv = cos*Lorentz_deriv_re + sin*Lorentz_deriv_im

    fig_scale = 1
    fontsize = 12

    plt.figure(figsize=(fig_scale*6.4,fig_scale*4.8))
    plt.xticks(fontsize=fontsize-2)
    plt.yticks(fontsize=fontsize-2)

    plt.xlim(-20.5,20.5)
    plt.ylim(-1.08,1.08)
    plt.yticks([-1,0,1])

    plt.plot(freq_axis,epr_J[:,1],'-',label="EPR spectrum")
    plt.plot(freq_axis,scale*Lorentz_deriv,'--',label="Lorentzian fit")

    plt.axvline(x=0,linestyle=':',color='gray')
    plt.axhline(y=0,linestyle=':',color='gray')

    plt.xlabel('Offset frequency ($\Omega/2\\pi$)  [MHz]',fontsize=fontsize)
    plt.ylabel('Intensity  [a.u.]',fontsize=fontsize)

    plt.legend(loc = 'upper right',fontsize=fontsize-2, facecolor='white',framealpha=1,frameon=False)

    plt.savefig("pdfs/p2_cw-epr.pdf", bbox_inches='tight')



##########
# for the fit to the experimental DNP spectrum
from scipy.optimize import minimize

def plot_dnp(dnp_enh,B1,T2e,T1e,R1n,b_nm,N_molar,tau,model,fit):

    delta2 = funs.delta2(b_nm, N_molar)

    fig_scale = 1
    fontsize = 14

    plt.rcParams['text.usetex'] = True

    freq_axis = funs.field_to_frequency(dnp_enh[:,0],'dnp')
    y_true = dnp_enh[:,1]

    #calculation
    freqI = 400
    omega_I = 2*np.pi* freqI

    R1n_del2 = R1n/delta2

    r1, r2, r3 = T1e/T2e, T2e/tau, 1

    
    if fit:
        def objective(x):
            r1, r2, r3 = x
            enh_calc = funs.calculate_enhancement(model,freq_axis,omega_I,B1,T2e,R1n_del2,r1,r2,r3)
            error = funs.RMSDerror(y_true,enh_calc)
            #print("error:",error)
            return error


        bounds = [(1,20),(1,100),(0.1,20)]
        initial = np.array([r1,r2,r3])
        res = minimize(objective,initial, method='SLSQP',bounds=bounds)
        print(res)
        r1, r2, r3 = res.x
    

    fit_tau = T2e/r2 # us
    fit_b = b_nm*np.power(1/r3,1/3) # nm
    print(f"T1:{r1*T2e:.2f} us, tau:{1e3*fit_tau:.2f} ns, b:{fit_b:.3} nm")
    print(f"B1*sqrt(r1):{B1*np.sqrt(r1):.3f}")
    print(f"D: {fit_b**2/fit_tau:.1f} nm2/us")

    enh_calc = funs.calculate_enhancement(model,freq_axis,omega_I,B1,T2e,R1n_del2,r1,r2,r3)
    #error = funs.RMSDerror(y_true,enh_calc)
    #print("error:",error)

    nrows, ncols = 2, 1

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', gridspec_kw={'height_ratios': [4, 3]},
                           figsize=(fig_scale*(ncols*5.+0.4),fig_scale*(nrows*2.4+0.4)))
    plt.subplots_adjust(hspace=0.05)


    for j in range(nrows):

        ax[j].axvline(x=-freqI,linestyle='--',color='gray')
        ax[j].axvline(x=0,linestyle=':',color='gray')
        ax[j].axvline(x=+freqI,linestyle='--',color='gray')
        ax[j].axhline(y=0,linestyle=':',color='gray')

        if j==0:
            ax[j].plot(freq_axis,y_true,'or',markerfacecolor='none',label="experiment")
            ax[j].plot(freq_axis,enh_calc,'-k',markerfacecolor='none',label=f"{model} (fit)")
            ax[j].plot(freq_axis,y_true-enh_calc,'-g',label="difference")
        
            ax[j].set_ylabel('DNP enhancement',fontsize=fontsize)
            ax[j].legend(loc = 'upper left',fontsize=fontsize-2, facecolor='white',framealpha=1,frameon=True)

        else:
            R2 = 1/T2e
            R1 = R2/r1
            tau = T2e/r2
            delta = np.sqrt(delta2*r3)
            #scan offset frequency
            pntN = 100
            Omega = 2*np.pi * np.linspace(freq_axis[0],freq_axis[-1],pntN) #M rad/s
            #Bloch
            omega1 = funs.omega1_from_B1(B1)
            fx, fy, p = funs.bloch_ss(R1,R2,omega1,Omega)
            sz =  R1 * p
            sx =  omega1 * fx * sz

            Tii_f,Tii0_f,vp_p, Tix_f, Tiy_f = funs.eig_numerical(R1,R2,omega1,omega_I,Omega,model,tau)
            Tiz_p = omega1*(fx*Tix_f - fy*Tiy_f)
            Tixp_p = Tiz_p/omega1/fx

            sx_norm = sx/np.max(sx)
            Tixp_norm = Tixp_p/np.max(Tixp_p)
            prod = sx*Tixp_p
            prod_norm = prod/np.max(prod[int(pntN/4):3*int(pntN/4)])

            ax[j].plot(Omega/(2*np.pi),sx_norm,'-',color='#1f77b4',label=r"$s_x^{\rm ss}$")#,color='#ff7f0e'
            
            ax[j].plot(Omega/(2*np.pi),Tixp_norm,color='#ff7f0e',label=r"$T_x$")
            ax[j].plot(Omega/(2*np.pi),prod_norm,'-.k',label=r"$s_x^{\rm ss} T_x$")
            ax[j].set_xlabel('Offset frequency ($\Omega/2\\pi$)  [MHz]',labelpad=10,fontsize=fontsize)

            ax[j].set_yticks([])
            ax[j].set_ylabel('a.u.',fontsize=fontsize,labelpad=28)

            ax[j].legend(loc = 'lower right',fontsize=fontsize-2, facecolor='white',framealpha=1,frameon=False,ncol=3)

    if model == 'ffhs':
        print("Paper 2, Figure 7\n")
        plt.savefig(f"pdfs/p2_fig7_B{B1}.pdf", bbox_inches='tight')
    elif model == 'exp':
        print("Paper 2, Figure 8\n")
        plt.savefig(f"pdfs/p2_fig8_B{B1}.pdf", bbox_inches='tight')
    
    return True 

