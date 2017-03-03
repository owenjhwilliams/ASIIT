'''
This is a set of functions for dealing with PIV data
Most functions are for 2D data. 
'''

#Import matlab PIV data
def importMatlabPIVdata(path,mats,structs):
    '''
    Plot 2D scalar fields
    
    Inputs: 
    path - Path to the matlab file
    mats - names of the matlab matrices to be imported.
    structs - names of the structures to be imported as dictionaries
    
    Output:
    Each of the matrices and structures in the order you specified in the input
    '''
    
    import numpy as np
    import h5py
    import pandas

    f = h5py.File(path)
    
    print(list(f.keys()))
    
    Temp = np.asarray(['Swirl'])
    
    ret = []

    for i in mats:
        Temp = np.asarray(f[i])
        if Temp.ndim == 2:
            Temp = np.transpose(Temp,(1,0))
        else:
            Temp = np.transpose(Temp,(2,1,0))
        ret.append(Temp)
        del Temp

    for i in structs:
        TempS = {k : f[i][k].value[0]       #Generate a dictionary linking all values in cond with their names
             for k in f[i].keys()}
        ret.append(TempS)
        del TempS
        
    return ret

#Plot 2D scalar field 
def plotScalarField(S,X=None,Y=None,bound=None,saveFolder=None):
    '''
    Plot 2D scalar fields
    
    Inputs: 
    X - 2D array with columns constant
    Y - 2D array with rows constant
    S - the 2D field to be plotted. Must be the same size and shape as X,Y
    bound = the symetric bound of the colorbar (-bound, bound), detault is max of abs(S)/5
    
    Output:
    Displays a plot of the scalar field
    '''
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    if bound is None:
        bound = np.round(np.max(np.absolute(S))/5)
        
    f = plt.figure(figsize = [8,3])
    if X is None:
        plt.pcolor(S, cmap='RdBu_r');
        plt.axis('scaled')
        plt.xlim([0, S.shape[1]])
        plt.ylim([0, S.shape[0]])
    else:
        plt.pcolor(X,Y,S, cmap='RdBu_r');
        plt.axis('scaled')
        plt.xlim([X.min(), X.max()])
        plt.ylim([Y.min(), Y.max()])
        
    plt.clim([-1*bound, bound])
    plt.colorbar()
    
    if saveFolder is not None:
        f.savefig(saveFolder, transparent=True, bbox_inches='tight', pad_inches=0)
    
    return [f, plt.gca()]


#Gets locations of distinct scalar blobs in each frame that are bigger than a certain threshold (in number of vectors)
def findBlobsSlow(S,Thresh=None):
    '''
    Note!!!! This code does not work with 3D data yet. 
    Finds distinct blobs of a scalar that are bigger than a certain size (Thresh)
    
    Inputs:
    S - sets of 2D scalar fields that have already been thresholded (0s or 1s)
    Thresh - Number of vectors that must be contained in a blob. If not defined, then no threshold filter will be used
    
    Outputs:
    cent
    labelled_array
    num_features
    features_per_frame
    
    '''
    import numpy as np
    from scipy.ndimage.measurements import label,find_objects,center_of_mass
    import copy
    
    uSize = S.shape
    
    labeled_array, num_features = label(S)
    print('There are ', num_features, ' features initially identified')
    
    if Thresh is not None:
        loc = find_objects(labeled_array)
        labeled_array_init = copy.copy(labeled_array)
        labeled_array[:] = 0;
        num_features_init = copy.copy(num_features)
        num_features = 0;
        for i in range(num_features_init):
            #print(np.max(labeled_array_init[loc[i]]),labeled_array_init[loc[i]],np.count_nonzero(labeled_array_init[loc[i]]))
            #print(labeled_array_init[loc[i]])
            #print(np.max(labeled_array_init[loc[i]]),np.count_nonzero(labeled_array_init[loc[i]]))
            if np.count_nonzero(labeled_array_init[loc[i]])>Thresh:
                #print('yes')
                num_features += 1;
                labeled_array[labeled_array_init==i+1] = num_features

        print('A total of ', num_features, ' are larger than the threshold size')
    
    features_per_frame = np.zeros(uSize[2],dtype=int);
    cent = [];
    for i in range(uSize[2]):
        features_per_frame[i] = len(np.unique(labeled_array[:,:,i])[1:])
        cent.append(center_of_mass(S[:,:,i],labeled_array[:,:,i],np.unique(labeled_array[:,:,i])[1:]))

    return [num_features, features_per_frame, labeled_array, cent]

#Gets locations of distinct scalar blobs in each frame that are bigger than a certain threshold (in number of vectors)
def findBlobs(S,Thresh=None,EdgeBound=None):
    '''
    Finds distinct blobs of a scalar that are bigger than a certain size (Thresh)
    Now new and improved! and much faster!
    
    Inputs:
    S - sets of 2D scalar fields that have already been thresholded (0s or 1s). The third dimension denotes the frame
    Thresh - Number of vectors that must be contained in a blob. If not defined, then no threshold filter will be used
    EdgeBound - Crops all blobs that are too close to the edge of the domain. No crop if left as none. 
    
    Outputs:
    cent - 
    labelled_array - The labeled array of blobs (in format of ndimage.measurements.label function). This is all the labels including the ones that might be too close to edge of domain
    num_features - Total number of features accross datasets
    features_per frame - Number of features identified in each frame
    
    '''
    import numpy as np
    from scipy.ndimage.measurements import label,find_objects,center_of_mass
    
    uSize = S.shape
    
    if S.ndim == 3:
        str_3D=np.array([[[0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]],

       [[0, 1, 0],
        [1, 1, 1],
        [0, 1, 0]],

       [[0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]]], dtype='uint8')
        #str_3D = str_3D.transpose(2,1,0)

        labeled_array, num_features = label(S.transpose(2,1,0),str_3D)
        labeled_array = labeled_array.transpose(2,1,0)
    else:
        labeled_array, num_features = label(S)
            
    #print(np.unique(labeled_array))
    #print(np.unique(labeled_array[:,:,0]))
    #print(labeled_array.shape)
    
    print('There are ', num_features, ' features identified')
    
    if Thresh is not None:
        loc = find_objects(labeled_array)
        labeled_array_out = labeled_array.copy()
        
        counts = np.bincount(labeled_array.ravel())
        
        ind = np.where(counts>Thresh)[0][1:]
        mask = np.in1d(labeled_array.ravel(), ind).reshape(labeled_array.shape)
        labeled_array_out[~mask] = 0
        
        [_, labeled_array_out] = np.unique(labeled_array_out,return_inverse=True)
        labeled_array_out = labeled_array_out.reshape(labeled_array.shape)
        
        num_features_out = len(ind)
        
        print('A total of ', num_features_out, ' are larger than the threshold size')
    else:
        labeled_array_out = labeled_array
        num_features_out = num_features
        
    #a = labeled_array_out[:,:,535]
    #b = a.copy()
    #c = list(range(uSize[1]))
    #b[:,:] = list(range(uSize[1])))
    #for i in np.unique(a):
        #print(np.unique(a)[15])
        #print(a[a==i])
            
    features_per_frame = np.zeros(uSize[2],dtype=int);
    cent = [];
    for i in range(uSize[2]):
        features_per_frame[i] = len(np.unique(labeled_array_out[:,:,i])[1:])
        cent.append(center_of_mass(S[:,:,i],labeled_array_out[:,:,i],np.unique(labeled_array_out[:,:,i])[1:]))
        
    #Round center locations to nearest index
    for i in range(len(cent)):
        for j in range(len(cent[i])):
            cent[i][j] = (int(round(cent[i][j][0])), int(round(cent[i][j][1])))
    
    #Remove all centers too close to edge of domain
    if EdgeBound is not None:
        newCent = []
        for i in range(len(cent)):
            newCent.append([])
            features_per_frame[i] = 0
            for j in range(len(cent[i])):
                if (cent[i][j][0]>EdgeBound-1 and cent[i][j][0]<uSize[0]-EdgeBound) and (cent[i][j][1]>EdgeBound-1 and cent[i][j][1] <uSize[1]-EdgeBound):
                    newCent[i].append(cent[i][j])
                    features_per_frame[i]+=1
        num_features_out = sum(features_per_frame)
        cent = newCent
        
        print('Of these', num_features_out, ' are far enough away from edge of domain')

    return [num_features_out, features_per_frame, labeled_array_out, cent]

# given the centers of a blobs, this function disects the vector field into a number of thumbnails of size frame x frame
def getThumbnails2D(Uf,Vf,Sf,cent,BoxSize):
    '''
    Given the centers of a blobs, this function disects the vector field into a number of thumbnails of size frame x frame. 
    All vectors inside frame but outside domain are padded with nans
    
    Inputs:
    
    Outputs:
    
    '''
    import numpy as np
    
    uSize = Uf.shape
    
    #find out how many features there are 
    #Round all centroids to integers
    num_features = 0
    for i in range(len(cent)):
        for j in range(len(cent[i])):
            #print(i, j)
            cent[i][j] = (int(round(cent[i][j][0])), int(round(cent[i][j][1])))
            num_features += 1
    
    #initialize thumbnail matrices
    Ut = np.zeros([2*BoxSize+1,2*BoxSize+1,num_features])    
    Ut[:] = np.NAN
    Vt = Ut.copy()
    St = Ut.copy()
    #print(Ut.shape)
    
    #pad out velocity fields so that there are NaNs around in all directions
    Uf2 = np.zeros([uSize[0]+2*BoxSize,uSize[1]+2*BoxSize,uSize[2]])    
    Uf2[:] = np.NAN
    Vf2 = Uf2.copy()
    Sf2 = Uf2.copy()

    Uf2[BoxSize:-1*BoxSize,BoxSize:-1*BoxSize,:] = Uf.copy()
    Vf2[BoxSize:-1*BoxSize,BoxSize:-1*BoxSize,:] = Vf.copy()
    Sf2[BoxSize:-1*BoxSize,BoxSize:-1*BoxSize,:] = Sf.copy()
    
    #print(Uf.shape)
    #print(Uf2.shape)
            
    #[f, ax] = PIVutils.plotScalarField(Sf2[:,:,0],bound=10)

    #for i in range(features_per_frame[1]):
        #plt.plot(cent[0][i][1]+BoxSize,cent[0][i][0]+BoxSize,'oy',markersize=4,markeredgecolor=None)
        
    #Now get the thumbnails
    thumb = 0
    for i in range(len(cent)):
        for j in range(len(cent[i])):
            Ut[:,:,thumb] = Uf2[cent[i][j][0]:cent[i][j][0]+2*BoxSize+1,cent[i][j][1]:cent[i][j][1]+2*BoxSize+1,i]  
            Vt[:,:,thumb] = Vf2[cent[i][j][0]:cent[i][j][0]+2*BoxSize+1,cent[i][j][1]:cent[i][j][1]+2*BoxSize+1,i] 
            St[:,:,thumb] = Sf2[cent[i][j][0]:cent[i][j][0]+2*BoxSize+1,cent[i][j][1]:cent[i][j][1]+2*BoxSize+1,i]  
            thumb+=1
            
    return [Ut, Vt, St]

# given the centers of a blobs, this function disects the vector field into a number of thumbnails of size frame x frame
def getThumbnails2D_noEdge(Uf,Vf,Sf,cent,BoxSize):
    '''
    Given the centers of a blobs, this function disects the vector field into a number of thumbnails of size frame x frame. 
    Domain is no longer padded with nans. A thumbnail is only taken if the center of the vortex is far enough from the frame    edges. 
    
    Inputs:
    
    Outputs:
    
    '''
    import numpy as np
    
    uSize = Uf.shape
    
    #find out how many features there are 
    #Round all centroids to integers
    num_features = 0
    for i in range(len(cent)):
        for j in range(len(cent[i])):
            #print(i, j)
            cent[i][j] = (int(round(cent[i][j][0])), int(round(cent[i][j][1])))
            num_features += 1
    
    #initialize thumbnail matrices
    Ut = np.zeros([2*BoxSize+1,2*BoxSize+1,num_features])    
    Ut[:] = np.NAN
    Vt = Ut.copy()
    St = Ut.copy()
    #print(Ut.shape)
    
    #Now get the thumbnails
    thumb = 0
    for i in range(len(cent)):
        for j in range(len(cent[i])):
            if (cent[i][j][0]>BoxSize-1 and cent[i][j][0]<uSize[0]-BoxSize) and (cent[i][j][1]>BoxSize-1 and cent[i][j][1] <uSize[1]-BoxSize):
                Ut[:,:,thumb] = Uf[cent[i][j][0]-BoxSize:cent[i][j][0]+BoxSize+1,cent[i][j][1]-BoxSize:cent[i][j][1]+BoxSize+1,i]
                Vt[:,:,thumb] = Vf[cent[i][j][0]-BoxSize:cent[i][j][0]+BoxSize+1,cent[i][j][1]-BoxSize:cent[i][j][1]+BoxSize+1,i]
                St[:,:,thumb] = Sf[cent[i][j][0]-BoxSize:cent[i][j][0]+BoxSize+1,cent[i][j][1]-BoxSize:cent[i][j][1]+BoxSize+1,i]  
                thumb+=1
            
    num_features = thumb-1
    return [Ut, Vt, St, num_features]

def getRandomThumbnails2D(Uf,Vf,Sf,numSamp,BoxSize):
    import numpy as np
    
    uSize = Uf.shape
    
    Pos = np.random.rand(3,numSamp)
    Pos[0] = Pos[0]*(uSize[0]-2*BoxSize-1)+BoxSize
    Pos[1] = Pos[1]*(uSize[1]-2*BoxSize-1)+BoxSize
    Pos[2] = Pos[2]* (uSize[2]-1)
    Pos = Pos.round().astype(int)
    
    #print(Pos[:,0])
    #print(np.ndarray.min(Pos[0]))
    #print(np.ndarray.max(Pos[0]))
    #print(np.ndarray.min(Pos[1]))
    #print(np.ndarray.max(Pos[1]))
    #print(np.ndarray.min(Pos[2]))
    #print(np.ndarray.max(Pos[2]))
    #print(max(Pos))
    
    #initialize thumbnail matrices
    Ut = np.zeros([2*BoxSize+1,2*BoxSize+1,numSamp])    
    Ut[:] = np.NAN
    Vt = Ut.copy()
    St = Ut.copy()
    #print(Ut.shape)
    
    #pad out velocity fields so that there are NaNs around in all directions
    Uf2 = np.zeros([uSize[0]+2*BoxSize,uSize[1]+2*BoxSize,uSize[2]])    
    Uf2[:] = np.NAN
    Vf2 = Uf2.copy()
    Sf2 = Uf2.copy()

    Uf2[BoxSize:-1*BoxSize,BoxSize:-1*BoxSize,:] = Uf.copy()
    Vf2[BoxSize:-1*BoxSize,BoxSize:-1*BoxSize,:] = Vf.copy()
    Sf2[BoxSize:-1*BoxSize,BoxSize:-1*BoxSize,:] = Sf.copy()
    
    #Now get the thumbnails
    thumb = 0
    for i in range(numSamp):
        #print(i)
        Ut[:,:,thumb] = Uf2[Pos[0,i]:Pos[0,i]+2*BoxSize+1,Pos[1,i]:Pos[1,i]+2*BoxSize+1,Pos[2,i]]  
        Vt[:,:,thumb] = Vf2[Pos[0,i]:Pos[0,i]+2*BoxSize+1,Pos[1,i]:Pos[1,i]+2*BoxSize+1,Pos[2,i]] 
        St[:,:,thumb] = Sf2[Pos[0,i]:Pos[0,i]+2*BoxSize+1,Pos[1,i]:Pos[1,i]+2*BoxSize+1,Pos[2,i]]  
        thumb+=1
            
    return [Ut, Vt, St]

def genHairpinField(BoxSize,Circ,r,rs,Ts,Rot,StagStren,Gvort,Gstag,Conv,x=None,y=None):
    '''
    Generates a theoretical hairpin vortex velocity field given a number of parameters. Returns U and V velocity fields
    
    Inputs:
    BoxSize - 2*Boxsize+1 is the number of vectors per side of box. 
    Circ - Circulation strength of vortex 
    r - diameter of vortex solid body rotation (constant vector magnitude outside core)
    rs, Ts - polar coordinate location of stagnation point
    Rot - Rotation of stagnation point shear layer
    StagStren - Vector magnitude of stagnation point velocity field (constant magnitude)
    Gvort - Width of blending gaussian for vortex
    Gstag - Width of blending gaussian for stagnation point
    Conv - Convective velocity of this vortex relative to the local mean
    
    Outputs:
    U - Streamwise velocity field
    V - Wall-normal velocity field
    '''
    
    import numpy as np
    import math
    
    #print((x.shape[0]-1)/2)
    
    if x is None:
        X, Y = np.meshgrid(np.arange(-1*BoxSize, BoxSize+1), np.arange(-1*BoxSize, BoxSize+1))
    else:
        assert BoxSize==(x.shape[0]-1)/2, 'The BoxSize does not match the length of the x vector. Thats not right...'
        assert BoxSize==(y.shape[0]-1)/2, 'The BoxSize does not match the length of the y vector. Thats not right...'
        X, Y = np.meshgrid(x, y)

    U = np.zeros([2*BoxSize+1,2*BoxSize+1])
    V = U.copy()
    R = np.hypot(X, Y)
    T = np.arctan2(Y,X)

    #Vortex
    Ut = Circ*R/(2*np.pi*r**2)
    #Ut[R>=r] = Circ/(2*np.pi*R[R>=r])
    Ut[R>=r] = Circ/(2*np.pi*r)      #make velocities constant outside core

    #Now convert back to cartesian velocities
    Uvort = Ut*np.sin(T)
    Vvort = -1*Ut*np.cos(T)

    #Create stagnation point flow
    Rot = Rot*np.pi/180*2
    Ts = Ts*np.pi/180
    xs = rs*np.cos(Ts)            #shift in stagnation point in x
    ys = rs*np.sin(Ts)            #shift in stagnation point in y
    #StagStren = 2;

    Xs = X-xs
    Ys = Y-ys;
    Ts = np.arctan2(Ys,Xs)

    M = np.hypot(Xs, Ys)
    U = M*np.cos(Ts-Rot)
    V = -1*M*np.sin(Ts-Rot)
    M = np.hypot(U, V)
    Ustag = U/M*StagStren
    Vstag = V/M*StagStren

    Ustag[np.isnan(U)] = 0
    Vstag[np.isnan(V)] = 0
    
    #Combine fields 
    Gvort_x = Gvort      #Radius of gaussian weighting function for vortex field
    Gvort_y = Gvort_x
    #Gstag = 5      #Radius of gaussian weighting function for stagnation point field

    Wvort = np.exp(-((X**2)/(2*Gvort_x)+(Y**2)/(2*Gvort_y)))
    Wvort_inv = -1*Wvort+1                #invert the weightings so that only vortex appears at vortex location                       

    Rstag = np.hypot(Xs, Ys)
    Wstag = np.exp(-((X-xs)**2/(2*Gstag)+(Y-ys)**2/(2*Gstag)))
    Wstag_inv = -1*Wstag+1

    #U = (Wvort*Uvort+Wstag*Ustag)/(Wvort+Wstag)
    #V = (Wvort*Vvort+Wstag*Vstag)/(Wvort+Wstag)

    U = (Wvort*Wstag_inv*Uvort+Wstag*Wvort_inv*Ustag)/(Wvort*Wstag_inv+Wstag*Wvort_inv)+Conv
    V = (Wvort*Wstag_inv*Vvort+Wstag*Wvort_inv*Vstag)/(Wvort*Wstag_inv+Wstag*Wvort_inv)

    #U = (Wstag_inv*Uvort+Wvort_inv*Ustag)/(Wvort_inv+Wstag_inv)
    #V = (Wstag_inv*Vvort+Wvort_inv*Vstag)/(Wvort_inv+Wstag_inv)
    
    return [U, V]