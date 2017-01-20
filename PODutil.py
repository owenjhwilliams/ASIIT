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
