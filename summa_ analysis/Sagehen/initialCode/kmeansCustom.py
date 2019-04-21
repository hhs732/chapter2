import matplotlib.pyplot as plt
import numpy as np

class K_Means:
    def __init__(self, numOfClusters=2, init_centroids=None):
        self.numOfClusters = numOfClusters
        self.centroids={}        
        for i in range(self.numOfClusters):
            self.centroids[i] = init_centroids[i]

    def fit(self,data):
        self.classifications = {}

        for i in range(self.numOfClusters):
            self.classifications[i] = []

        for featureset in data:
            distances = [np.linalg.norm(featureset-self.centroids[centroid]) for centroid in self.centroids]
            classification = distances.index(min(distances))
            self.classifications[classification].append(featureset)


    def predict(self,data):
        distances = [np.linalg.norm(data-self.centroids[centroid]) for centroid in self.centroids]
        classification = distances.index(min(distances))
        return classification


#%%
   
# this is your data     
X = np.array([[1, 2],
              [1.5, 1.8],
              [5, 8 ],
              [8, 8],
              [1, 0.6],
              [9,11],
              [1,3],
              [8,9],
              [0,3],
              [5,4],
              [6,4],])


# this is given centroids
centroids=np.array([[1,1],[2,2],[3,3]])    

# instantiate a class
clf = K_Means(numOfClusters=3,init_centroids=centroids)

# fit kmean class to data
clf.fit(X)

# get classification 
print clf.classifications
print ("-"*50)

# get predictions
for x in X:
    y=clf.predict(x)
    print("%s is close to  %s" %(x,centroids[y]))