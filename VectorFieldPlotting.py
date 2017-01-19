'''
This is a set of plotting functions to plot vector fields from PIV data
Most functions are for 2D data
'''

#Plot 2D scalar field 
def plotScalarField(X,Y,S,bound):
    '''
    Plot 2D scalar fields
    
    Inputs: 
    X - 2D array with columns constant
    Y - 2D array with rows constant
    S - the 2D field to be plotted. Must be the same size and shape as X,Y
    bound = the symetric bound of the colorbar (-bound, bound)
    
    Output:
    Plots the modes corresponding to plotModes
    '''
    
    plt.figure(figsize = [8,3])
    plt.pcolor(X,Y,S, cmap='RdBu_r');
    plt.clim([-1*bound, bound])
    plt.axis('scaled')
    plt.xlim([X.min(), X.max()])
    plt.ylim([Y.min(), Y.max()])
    plt.colorbar()