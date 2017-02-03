''' 
This is a set of utility funcitons useful for analysing POD data. Plotting and data reorganization functions
'''

#Plot 2D POD modes 
def plotPODmodes2D(X,Y,Umodes,Vmodes,plotModes):
    '''
    Plot 2D POD modes
    
    Inputs: 
    X - 2D array with columns constant
    Y - 2D array with rows constant
    Umodes - 3D array with Umodes.shape[2] = the total number of modes plotted
    Vmodes - 3D array with Umodes.shape[2] = the total number of modes plotted
    plotModes = a list of modes to be plotted
    
    Output:
    Plots the modes corresponding to plotModes
    '''

    import matplotlib.pyplot as plt
    
    #assert plotModes.max()<=Umodes.shape[2], 'You asked for more modes than were calculated'
    assert Umodes.shape[2]==Vmodes.shape[2], 'There are different numbers of U and V modes. Thats not right...'
    
    for i in plotModes:
        f, ax = plt.subplots(1,2)
        f.set_figwidth(18)
        im1 = ax[0].pcolor(X,Y,Umodes[:,:,i],cmap='RdBu_r');
        im2 = ax[1].pcolor(X,Y,Vmodes[:,:,i],cmap='RdBu_r');

        ax[0].set_title('U - #' + str(i+1))
        ax[0].set_aspect('equal')
        ax[0].set_xlim([X.min(),X.max()])
        ax[0].set_ylabel('y(m)')
        ax[0].set_xlabel('x(m)')
        
        ax[1].set_title('V - #' + str(i+1))
        ax[1].set_aspect('equal')
        ax[1].set_xlim([X.min(),X.max()])
        ax[1].set_ylabel('y(m)')
        ax[1].set_xlabel('x(m)')

        cbar1 = f.colorbar(im1,ax=ax[0])
        im1.set_clim(-1*max(map(abs,cbar1.get_clim())), max(map(abs,cbar1.get_clim()))) 
        cbar2 = f.colorbar(im2,ax=ax[1])
        im2.set_clim(-1*max(map(abs,cbar2.get_clim())), max(map(abs,cbar2.get_clim()))) 
        
        del im1,im2,cbar1,cbar2

#Plot 3D POD modes 
def plotPODmodes3D(X,Y,Umodes,Vmodes,Wmodes,plotModes):
    '''
    Plot 2D POD modes
    
    Inputs: 
    X - 2D array with columns constant
    Y - 2D array with rows constant
    Umodes - 3D array with Umodes.shape[2] = the total number of modes plotted
    Vmodes - 3D array with Umodes.shape[2] = the total number of modes plotted
    Wmodes - 3D array with Umodes.shape[2] = the total number of modes plotted
    plotModes = a list of modes to be plotted
    
    Output:
    Plots the modes corresponding to plotModes
    '''

    import matplotlib.pyplot as plt
    
    #assert plotModes.max()<=Umodes.shape[2], 'You asked for more modes than were calculated'
    assert Umodes.shape[2]==Vmodes.shape[2], 'There are different numbers of U and V modes. Thats not right...'
    assert Umodes.shape[2]==Wmodes.shape[2], 'There are different numbers of U and W modes. Thats not right...'
    
    for i in plotModes:
        f, ax = plt.subplots(1,3)
        f.set_figwidth(18)
        im1 = ax[0].pcolor(X,Y,Umodes[:,:,i],cmap='RdBu_r');
        im2 = ax[1].pcolor(X,Y,Vmodes[:,:,i],cmap='RdBu_r');
        im3 = ax[2].pcolor(X,Y,Wmodes[:,:,i],cmap='RdBu_r');

        ax[0].set_title('U - #' + str(i+1))
        ax[0].set_aspect('equal')
        ax[0].set_xlim([X.min(),X.max()])
        ax[0].set_ylabel('y(m)')
        ax[0].set_xlabel('x(m)')
        
        ax[1].set_title('V - #' + str(i+1))
        ax[1].set_aspect('equal')
        ax[1].set_xlim([X.min(),X.max()])
        ax[1].set_ylabel('y(m)')
        ax[1].set_xlabel('x(m)')
        
        ax[2].set_title('W - #' + str(i+1))
        ax[2].set_aspect('equal')
        ax[2].set_xlim([X.min(),X.max()])
        ax[2].set_ylabel('y(m)')
        ax[2].set_xlabel('x(m)')

        cbar1 = f.colorbar(im1,ax=ax[0])
        im1.set_clim(-1*max(map(abs,cbar1.get_clim())), max(map(abs,cbar1.get_clim()))) 
        cbar2 = f.colorbar(im2,ax=ax[1])
        im2.set_clim(-1*max(map(abs,cbar2.get_clim())), max(map(abs,cbar2.get_clim()))) 
        cbar3 = f.colorbar(im3,ax=ax[2])
        im3.set_clim(-1*max(map(abs,cbar3.get_clim())), max(map(abs,cbar3.get_clim()))) 
        
        del im1,im2,im3,cbar1,cbar2,cbar3

#Reorganize modes matrix so that the modes can be easily plotted
def reconstructPODmodes(modes,uSize,num_modes,numC):
    '''
    Reconstruct the mode shapes for three component single plane data
    
    Inputs: 
    modes - outout from mr.compute_POD_matrices_snaps_method
    uSize - size of original velocity dataset
    num_modes - number of modes calculated by mr.compute_POD_matrices_snaps_method
    numC - number of velocity components
    
    Output:
    Umodes, Vmodes and optionally Wmodes
    '''       
    import numpy as np
    
    #Rearrange mode data to get mode fields
    modeSize = modes.shape
    Umodes = modes[0:uSize[0]*uSize[1],:];
    Umodes2 = np.zeros((uSize[0],uSize[1],num_modes))
    
    if numC >= 2:
        Vmodes = modes[uSize[0]*uSize[1]:2*uSize[0]*uSize[1],:];
        Vmodes2 = np.zeros((uSize[0],uSize[1],num_modes))
    if numC >= 3:
        Wmodes = modes[2*uSize[0]*uSize[1]:modeSize[0]+1,:];
        Wmodes2 = np.zeros((uSize[0],uSize[1],num_modes))

    Umodes.shape
    

    for i in range(num_modes):
        #i=1
        Umodes2[:,:,i] = np.reshape(Umodes[:,i],(uSize[0],uSize[1]))
        if numC >=2:
            Vmodes2[:,:,i] = np.reshape(Vmodes[:,i],(uSize[0],uSize[1]))
        if numC >=3:
            Wmodes2[:,:,i] = np.reshape(Vmodes[:,i],(uSize[0],uSize[1]))        

        #Umodes.shape
        #uSize[0]*uSize[1]
    if numC == 1:
        return [Umodes2]
    elif numC == 2:
        return [Umodes2, Vmodes2]
    elif numC == 3:
        return [Umodes2, Vmodes2, Wmodes2]
        
#Plot heatmaps of POD coefficients
def plotPODcoeff(C,modes,num_bins,bound=None,logscale=None):
    '''
    Reconstruct the mode shapes for three component single plane data
    
    Inputs: 
    C - matrix of coefficients (mode number, coefficent for each frame) 
    modes - indices of modes to be plotted 
    num_bins - size of bins to be plotted. Passed to hexbin
    bound - the axis bound. If none taken to be max coefficient
    logscale - describing whether or not to do the heatmap in a log scale
    
    Output:
    plots a grid of hexbin plots for each mode
    '''       
    import numpy as np
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt
    
    if bound == None:
        bound = round(np.max(np.absolute(C)))
    
    xedges = np.linspace(-1*bound, bound, num=num_bins)
    yedges = xedges;
    bound = 0.5*(xedges[1]+xedges[2]);
    
    Z, xedges, yedges = np.histogram2d(C[0], C[1], bins=(xedges, yedges))
    xv, yv = np.meshgrid(0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1]))

    fig, axs = plt.subplots(ncols=len(modes)-1,nrows=len(modes)-1,figsize=(9, 12))
    fig.subplots_adjust(hspace=0.01, left=0.01, right=1)
    
    #print(axs.shape)

    for i in range(len(modes)-1):
        for j in range(len(modes)-1):
            ax = axs[i,j]
            if j>=i:
                Z, xedges, yedges = np.histogram2d(C[i],C[j+1], bins=(xedges, yedges))
                if logscale == None:
                    hb = ax.pcolor(xv, yv, Z, cmap='hsv')
                else:
                    hb = ax.pcolor(xv, yv, np.log(Z+1), cmap='hsv')
                    
                ax.plot([-1*bound, bound],[0, 0],'--k')
                ax.plot([0, 0],[-1*bound, bound],'--k')
                
                if i == 0:
                    ax.set_xlabel('C{0}'.format(j+2))
                    ax.xaxis.tick_top()
                    ax.xaxis.set_label_position("top")
                    ax.tick_params(axis='x', labelsize=7)
                else:
                    ax.set_xticklabels([])

                    
                if j == len(modes)-2:
                    ax.yaxis.tick_right()
                    ax.set_ylabel('C{0}'.format(i+1))
                    ax.yaxis.set_label_position("right")
                    ax.tick_params(axis='y', labelsize=7)
                else:
                    ax.set_yticklabels([])
                    
                ax.set_xlim(bound,-1*bound)
                ax.set_ylim(bound,-1*bound)
                
                ax.set_aspect("equal")
                ax.set_adjustable("box-forced")

                #fig = plt.figure(figsize = [8,3])
                #hb = ax.hexbin(C[0], C[1], gridsize=10, cmap='OrRd')
                #plt.axis([-1*bound, bound, -1*bound, bound])
                #plt.axis('scaled')
                #cb = fig.colorbar(hb, ax=ax)
                #cb.set_label('counts')
                #cb.set_label('log10(N)')
            else:
                ax.axis('off')
           
        
#Plot heatmaps of POD coefficients
def plotYposPODcoeff(ypos,C,modes,num_bins,bound=None,logscale=None):
    '''
    Reconstruct the mode shapes for three component single plane data
    
    Inputs: 
    ypos - wall-normal position of each thumbnail. [0 1]
    C - matrix of coefficients (mode number, coefficent for each frame) 
    modes - indices of modes to be plotted 
    num_bins - size of bins to be plotted. Passed to hexbin
    bound - the axis bound. If none taken to be max coefficient
    logscale - describing whether or not to do the heatmap in a log scale
    
    Output:
    plots a grid of hexbin plots for each mode
    '''       
    import numpy as np
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt
    
    if bound == None:
        bound = round(np.max(np.absolute(C)))
    
    #xedges = np.linspace(0, 1, num=num_bins);
    xedges = ypos
    xedges = np.concatenate([[xedges[0]-(xedges[1]-xedges[0])],xedges, [xedges[-1]+xedges[-1]-xedges[-2]]])
    #print(xedges)
    yedges = np.linspace(-1*bound, bound, num=num_bins);
    #bound = 0.5*(xedges[1]+xedges[2]);
    
    Z, xedges, yedges = np.histogram2d(C[0],C[1], bins=(xedges, yedges))
    xv, yv = np.meshgrid(0.5*(xedges[1:]+xedges[:-1]),0.5*(yedges[1:]+yedges[:-1]))
    
    fig, axs = plt.subplots(nrows=len(modes),figsize=(3, 3*len(modes)))
    fig.subplots_adjust(hspace=0.1, left=0.1, right=1)
    
    #print(axs.shape)

    for i in range(len(modes)):
        ax = axs[i]

        Z, xedges, yedges = np.histogram2d(C[0],C[i+1], bins=(xedges, yedges))
        Z = Z.T
        if logscale == None:
            hb = ax.pcolor(xv, yv, Z, cmap='hsv')
        else:
            hb = ax.pcolor(xv, yv, np.log(Z+1), cmap='hsv')

        ax.plot([-1*bound, bound],[0, 0],'--k')

        if i == len(modes)-1:
            ax.set_xlabel('y/delta')
            ax.tick_params(axis='x', labelsize=7)
        else:
            ax.set_xticklabels([])

        ax.set_ylabel('C{0}'.format(i+1))
        ax.tick_params(axis='y', labelsize=7)

        ax.set_xlim(0,max(ypos))
        ax.set_ylim(-1*bound,bound)


